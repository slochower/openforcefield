#!/usr/bin/env python3

# To do
# Add something about spitting script version into output
# If only running a single phase, then we need to make sure the proper files are in place from earlier
# Allow the production phase to start at an arbitrary window

import numpy as np
import subprocess as sp
import re as re
import datetime as dt
import os as os
import socket as socket
import sys as sys
import glob as glob

class simulation():

    def __init__(self):
        # System parameters
        self.system_name = 'full-HMR'
        self.minimize = True
        self.thermalize = False
        self.equilibrate = False
        self.produce = False
        self.resume = False
        # File parameters
        self.minimization_prefix   = 'mini'
        self.thermalization_prefix = 'therm'
        self.equilibration_prefix  = 'eqnpt'
        self.production_prefix     = 'traj'
        self.log_file = 'system-test'
        self.log = None                         # this is necessary before we open the file for writing later
        # MD engines
        self.minimization_engine   = '$AMBERHOME/bin/pmemd'
        self.thermalization_engine = 'mpiexec.hydra -np 2 $AMBERHOME/bin/pmemd.cuda.MPI'
        self.equilibration_engine  = 'mpiexec.hydra -np 2 $AMBERHOME/bin/pmemd.cuda.MPI'
        self.production_engine     = 'mpiexec.hydra -np 2 $AMBERHOME/bin/pmemd.cuda.MPI'
        # self.cuda_visible_devices  = '0,1'
        self.cuda_visible_devices  = None       # can be a single number like '0', a list like '0,1',
                                                # 'tscc' where it gets set automatically, or None
        self.cuda_visible_devices  = 'tscc'
        # Minimization parameters
        self.minimization_cycles  = 20000       # maximum minimization cycles
        self.minimization_switch  = 5           # switch to CG from SD after $ cycles
        # Thermalization parameters
        self.thermalization_dt    = 0.004       # initial temperature
        self.thermalization_t0    = 10.0        # initial temperature
        self.thermalization_t1    = 300.0       # final temperature
        self.thermalization_cut   = 9           # angstroms
        self.thermalization_steps = 250000      # number of steps (nstlim)
        # Equilibration parameters
        self.equilibration_dt     = 0.004       # femtoseconds
        self.equilibration_temp   = 300         # kelvin
        self.equilibration_steps  = 250000      # number of MD steps (nstlim) per window
        self.equilibration_out    = 500         # trajectory output frequency (steps)
        self.equilibration_time   = 10          # total nanoseconds
        # Production MD NPT parameters
        self.production_dt     = 0.004       # femtoseconds
        self.production_temp   = 300         # kelvin
        self.production_steps  = 2500000     # number of MD steps (nstlim) per window
        self.production_out    = 500         # trajectory output frequency (steps)
        self.production_time   = 1000        # nanoseconds
        self.production_water  = False       # include water in trajectory
        self.production_ntwprt = 3431        # include this many atoms in trajectory (if water is False, then this
                                             # should be the number of atoms to keep; if water is True, then this option
                                             # is not used )

    def simulate(self):
        self.start_logging()
        self.print_parameters()
        if self.minimize:
            self.minimization()
        if self.thermalize:
            self.thermalization()
        if self.equilibrate:
            self.equilibration()
        if self.produce:
            self.production()
        self.end_logging()

    def now(self):
        d_date = dt.datetime.now()
        return d_date.strftime("%Y-%m-%d %I:%M:%S %p")


    def minimization(self):
        print('# Minimization for {} steps'.format(self.minimization_cycles), file=self.log)
        # Check if output file already exists
        if os.path.isfile(self.minimization_prefix + '.out'):
            with open(self.minimization_prefix + '.out', 'r') as f:
                for line in f:
                    if re.search("TIMING", line):
                        sys.exit('Minimization exists and it looks correct. Aborting.' )
        print('{:<25} {:<10}'.format('Initializing', self.now()), file=self.log)
        # Craft the input file
        string = '''
Minimizing.
 &cntrl
  imin = 1,         ! minimize
  ntx = 1,          ! coordinates will be read
  ntpr = 50,        ! number of steps to output human-readable energy
  maxcyc = {},      ! maximum number of cycles
  ncyc = {},        ! number of steps until switch from SD to CG
  irest = 0,        ! not a restart
  ntf = 1,          ! all interactions calculated
  ntc = 1,          ! shake is not performed
  ntb = 1,          ! constant volume (constant pressure n/a for minimization)
  igb = 0,          ! no GB or PB implicit solvent
  cut = 9.0,        ! nonbonded cutoff
 /

'''.format(self.minimization_cycles, self.minimization_switch)
        mini = open(self.minimization_prefix + '.in', 'w')
        mini.write(string)
        mini.close()
        # Form the command
        command = '{} -O -p {} -ref {} -c {} -i {} -o {} -r {} -inf /dev/null'.format(
            self.minimization_engine,
            self.system_name + '.prmtop',
            self.system_name + '.inpcrd',
            self.system_name + '.inpcrd',
            self.minimization_prefix + '.in',
            self.minimization_prefix + '.out',
            self.minimization_prefix + '.rst'
        )
        print('{:<25} {:<10}'.format('Running', self.now()), file=self.log)
        script = open(self.minimization_prefix + '.sh', 'w')
        script.write('#!/usr/bin/env bash\n')
        script.write('source $AMBERHOME/amber.sh\n')
        script.write(command)
        script.close()
        # Execute
        sp.call(['bash ./{}'.format(self.minimization_prefix + '.sh')], shell=True)
        # Check
        print('{:<25} {:<10}'.format('...checking', self.now()), file=self.log)
        try:
            sp.check_output(['grep -q TIMINGS {}'.format(self.minimization_prefix + '.out')], shell=True)
            print('{:<25} {:<10}'.format('...done', self.now()), file=self.log)
        except sp.CalledProcessError as e:
            print('{:<25} {:<10}'.format('...failed', self.now()), file=self.log)
            sys.exit('Uh oh: looks like minimization failed.')


    def thermalization(self):
        print('# NVT Thermalization for {} ns'.format(self.thermalization_steps * self.thermalization_dt / 1000), file=self.log)
        # Check if output file already exists
        if os.path.isfile(self.thermalization_prefix + '.out'):
            with open(self.thermalization_prefix + '.out', 'r') as f:
                for line in f:
                    if re.search("TIMING", line):
                        sys.exit('Thermalization exists and it looks correct. Aborting.' )
        print('{:<25} {:<10}'.format('Initializing', self.now()), file=self.log)
        # Craft the input file
        string = '''
Thermalizing, NVT.
 &cntrl
  imin = 0,         ! run without minimization
  ntx = 1,          ! coordinates but not velocities will be read
  irest = 0,        ! not a restart
  ntpr = 500,       ! number of steps to output human-readable energy
  ntwr = 500,       ! number of steps to output binary restart file
  ntwx = 0,         ! number of steps to output coordinates to trajectory
  ioutfm = 1,       ! netcdf
  ntf = 2,          ! bond interactions with H are omitted
  ntc = 2,          ! shake for hydrogen bonds
  cut = 9.0,        ! nonbonded cutoff
  ntt = 3,          ! langevin dynamics using gamma_ln (canonical ensemble)
  tempi = {}        ! temperature
  gamma_ln = 1.0,   ! collision frequency for ntt = 3
  ig = -1,          ! random number seed
  ntp = 0,          ! no pressue scaling
  nstlim = {},      ! number of md steps
  dt = {},          ! time step
  ntb = 1,          ! constant volume 
  igb = 0,          ! no GB or PB implicit solvent
  nmropt = 1,       ! restraints (below)
 /
 &wt type='TEMP0', istep1=0,    istep2={},  value1={}, value2={}, /
 &wt type = 'END', /
  
'''.format(self.thermalization_t0, self.thermalization_steps,
           self.thermalization_dt, self.thermalization_steps, self.thermalization_t0, self.thermalization_t1)
        therm = open(self.thermalization_prefix + '.in', 'w')
        therm.write(string)
        therm.close()
        # Form the command
        command = '{} -O -p {} -ref {} -c {} -i {} -o {} -r {} -x {} -inf {}'.format(
            self.thermalization_engine,
            self.system_name + '.prmtop',
            self.system_name + '.inpcrd',
            self.minimization_prefix + '.rst',
            self.thermalization_prefix + '.in',
            self.thermalization_prefix + '.out',
            self.thermalization_prefix + '.rst',
            self.thermalization_prefix + '.nc',
            self.thermalization_prefix + '.mdinfo',
        )
        print('{:<25} {:<10}'.format('Running', self.now()), file=self.log)
        script = open(self.thermalization_prefix + '.sh', 'w')
        script.write('#!/usr/bin/env bash\n')
        script.write('source $AMBERHOME/amber.sh\n')
        if self.cuda_visible_devices is not (None or 'tscc'):
            script.write('export CUDA_VISIBLE_DEVICES={}\n'.format(str(self.cuda_visible_devices)))
        script.write(command)
        script.close()
        # Execute
        sp.call(['bash ./{}'.format(self.thermalization_prefix + '.sh')], shell=True)
        # Check
        print('{:<25} {:<10}'.format('...checking', self.now()), file=self.log)
        try:
            sp.check_output(['grep -q TIMINGS {}'.format(self.thermalization_prefix + '.out')], shell=True)
            print('{:<25} {:<10}'.format('...done', self.now()), file=self.log)
        except sp.CalledProcessError as e:
            print('{:<25} {:<10}'.format('...failed', self.now()), file=self.log)
            sys.exit('Uh oh: looks like thermalization failed.')

    def equilibration(self):
        ns_per_file = self.equilibration_steps * self.equilibration_dt / 1000
        number_eq_files = int(self.equilibration_time / ns_per_file)
        print('# NPT Equilibration for {} ns in {} windows'.format(self.equilibration_time,
                                                                   number_eq_files), file=self.log)
        # Check if output file already exists
        if os.path.isfile(self.equilibration_prefix + '.out'):
            with open(self.equilibration_prefix + '.out', 'r') as f:
                for line in f:
                    if re.search("TIMING", line):
                        sys.exit('Equilibration exists and it looks correct. Aborting.' )
        print('{:<25} {:<10}'.format('Initializing', self.now()), file=self.log)
        # Craft the input file
        string = '''
Equilibrate {} ns, NPT.
 &cntrl
  imin = 0,            ! run without minimization
  ntx = 5,             ! coordinates and velocities will be read in
  irest = 1,           ! this is a restart
  ntpr = {},           ! number of steps to output human-readable energy
  ntwr = {},           ! number of steps to output binary restart file
  ntwx = {},           ! number of steps to output coordinates to trajectory
  ioutfm = 1,          ! netcdf
  ntf = 2,             ! bond interactions with H are omitted
  ntc = 2,             ! shake for hydrogen bonds
  cut = 9.0,           ! nonbonded cutoff 
  ntt = 3,             ! langevin dynamics using gamma_ln (canonical ensemble)
  temp0 = {},          ! temperature
  gamma_ln = 1.0,      ! collision frequency
  ig = -1,             ! random number seed
  nstlim = {},         ! number of md steps
  ntb = 2,             ! constant pressure 
  ntp = 1,             ! isotropic pressure scaling
  barostat = 2,        ! monte carlo barostat
  igb = 0,             ! no GB or PB implicit solvent
  dt = {},             ! time step
 / 
  
'''.format(ns_per_file, self.equilibration_out, self.equilibration_out, self.equilibration_out,
           self.equilibration_temp, self.equilibration_steps, self.equilibration_dt)
        eq = open(self.equilibration_prefix + '.in', 'w')
        eq.write(string)
        eq.close()
        # Now we need to write a loop to handle multiple windows...
        sp.call(['cp {} {}'.format(self.thermalization_prefix + '.rst',
                                  self.equilibration_prefix + '.0000.rst')], shell=True)
        for window in range(1, number_eq_files + 1):
            command = '{} -O -p {} -ref {} -c {} -i {} -o {} -r {} -x {} -inf /dev/null'.format(
                self.equilibration_engine,
                self.system_name + '.prmtop',
                self.system_name + '.inpcrd',
                self.equilibration_prefix + '.{:04d}.rst'.format(window-1),
                self.equilibration_prefix + '.in',
                self.equilibration_prefix + '.{:04d}.out'.format(window),
                self.equilibration_prefix + '.{:04d}.rst'.format(window),
                self.equilibration_prefix + '.{:04d}.nc'.format(window),
            )
            print('{:<25} {:<10}'.format('Window ' + str(window), self.now()), file=self.log)
            script = open(self.equilibration_prefix + '.{:04d}.sh'.format(window), 'w')
            script.write('#!/usr/bin/env bash\n')
            script.write('source $AMBERHOME/amber.sh\n')
            if self.cuda_visible_devices is not None:
                script.write('export CUDA_VISIBLE_DEVICES={}\n'.format(str(self.cuda_visible_devices)))
            script.write(command)
            script.close()
            # Execute
            sp.call(['bash ./{}.{:04d}.sh'.format(self.equilibration_prefix, window)], shell=True)
            # Check
            print('{:<25} {:<10}'.format('...checking', self.now()), file=self.log)
            try:
                sp.check_output(['grep -q TIMINGS {}.{:04d}.out'.format(self.equilibration_prefix, window)], shell=True)
                print('{:<25} {:<10}'.format('...done', self.now()), file=self.log)
            except sp.CalledProcessError as e:
                print('{:<25} {:<10}'.format('...failed', self.now()), file=self.log)
                sys.exit('Uh oh: looks like equilibration failed.')
            # Save the final restart
            sp.call(['cp {}.{:04d}.rst {}'.format(self.equilibration_prefix, window,
                                  self.production_prefix + '.0000.rst')], shell=True)

    def production(self):
        if self.resume:
            print('{:<25} {:<10}'.format('Resuming', self.now()), file=self.log)
            # Find the last trajectory file (whether completed or not)
            trajectories = sorted(glob.glob(self.production_prefix + '*.nc'))
            last_trajectory = trajectories[-1]
            # If we have more than one trajectory, we should start from the 2nd to last
            if len(trajectories) > 1:
                penultimate_trajectory = trajectories[-2]
                print('{:<25} {:<10}'.format('Checking', penultimate_trajectory), file=self.log)
                penultimate_window = penultimate_trajectory.rsplit(".", 1)[0]
                # Check the output from the penultimate trajectory is okay
                try:
                    grep = sp.check_output(['grep -q TIMINGS {}'.format(str(penultimate_window) + '.out')], shell=True)
                except sp.CalledProcessError as grep_exception:
                    print(grep_exception.returncode, grep_exception.output, file=self.log)
                    sys.exit('Penultimate trajectory {} did not complete successfully.'.format(
                        str(penultimate_window) + '.out'))
                # Make sure there is a restart file after the penultimate trajectory
                if not os.path.isfile(str(penultimate_window) + '.rst'):
                    sys.exit('Can\'t find penultimate trajectory restart file.')
                print('{:<25} {:<10}'.format('...done', self.now()), file=self.log)
                # Roll back to the penultimate trajectory file
                last_good_window = penultimate_window.rsplit(".", 1)[1]
                current_window = int(last_good_window) + 1
            elif len(trajectories) == 1:
                print('# Found a single good trajectory.', file=self.log)
                current_file = last_trajectory.rsplit(".", 1)[0]
                current_window = current_file.rsplit(".", 1)[1]
            else:
                print('# Error finding a good starting point.', file=self.log)

        else:
            current_window = 1

        ns_per_file = self.production_steps * self.production_dt / 1000
        number_production_files = int(self.production_time / ns_per_file)
        print('# NPT production for {} ns (total) in {} windows'.format(self.production_time,
                                                                number_production_files), file=self.log)
        print('# {} ns in each window'.format(self.production_dt / 1000 * self.production_steps), file=self.log)
        print('# {} frames in each window'.format(self.production_steps / self.production_out), file=self.log)
        print('# {} ps per frame in each window'.format(self.production_dt * self.production_out), file=self.log)
        # Check if output file already exists
        if os.path.isfile(self.production_prefix + '.out'):
            with open(self.production_prefix + '.out', 'r') as f:
                for line in f:
                    if re.search("TIMING", line):
                        sys.exit('production exists and it looks correct. Aborting.' )
        print('{:<25} {:<10}'.format('Initializing', self.now()), file=self.log)
        # Craft the input file
        if not self.production_water:
            string = '''
Production {} ns, NPT.
 &cntrl
  imin = 0,             ! run without minimization
  ntx = 5,              ! coordinates and velocities will be read in
  irest = 1,            ! this is a restart
  ntpr = {},            ! number of steps to output human-readable energy
  ntwr = {},            ! number of steps to output binary restart file
  ntwe = {},            ! number of steps to output binary energy and temperature
  ntxo = 2,             ! netcdf
  ntwx = {},            ! number of steps to output coordinates to trajectory
  ioutfm = 1,           ! netcdf
  ntf = 2,              ! bond interactions with H are omitted
  ntc = 2,              ! shake for hydrogen bonds
  cut = 9.0,            ! nonbonded cutoff
  ntt = 3,              ! langevin dynamics using gamma_ln (canonical ensemble)
  temp0 = {},           ! temperature
  gamma_ln = 1.0,       ! collision frequency
  ig = -1,              ! random number seed
  nstlim = {},          ! number of md steps
  dt = {},              ! time step
  ntb = 2,              ! constant pressure
  ntp = 1,              ! isotropic pressure scaling
  barostat = 2,         ! monte carlo barostat
  igb = 0,              ! no GB or PB implicit solvent
  ntwprt = {}           ! only include these atoms in trajectory
 /

'''.format(ns_per_file, self.production_out, self.production_out, int(self.production_out/2),
           self.production_out, self.production_temp, self.production_steps, self.production_dt,
           self.production_ntwprt)
        else:
            string = '''
            Production {} ns, NPT.
             &cntrl
              imin = 0,             ! run without minimization
              ntx = 5,              ! coordinates and velocities will be read in
              irest = 1,            ! this is a restart
              ntpr = {},            ! number of steps to output human-readable energy
              ntwr = {},            ! number of steps to output binary restart file
              ntwe = {},            ! number of steps to output binary energy and temperature
              ntxo = 2,             ! netcdf
              ntwx = {},            ! number of steps to output coordinates to trajectory
              ioutfm = 1,           ! netcdf
              ntf = 2,              ! bond interactions with H are omitted
              ntc = 2,              ! shake for hydrogen bonds
              cut = 9.0,            ! nonbonded cutoff
              ntt = 3,              ! langevin dynamics using gamma_ln (canonical ensemble)
              temp0 = {},           ! temperature
              gamma_ln = 1.0,       ! collision frequency
              ig = -1,              ! random number seed
              nstlim = {},          ! number of md steps
              dt = {},              ! time step
              ntb = 2,              ! constant pressure
              ntp = 1,              ! isotropic pressure scaling
              barostat = 2,         ! monte carlo barostat
              igb = 0,              ! no GB or PB implicit solvent
             /

            '''.format(ns_per_file, self.production_out, self.production_out, int(self.production_out / 2),
                       self.production_out, self.production_temp, self.production_steps, self.production_dt)

        prod = open(self.production_prefix + '.in', 'w')
        prod.write(string)
        prod.close()
        # Now we need to write a loop to handle multiple windows...
        for window in range(current_window, number_production_files + 1):
            command = '{} -O -p {} -ref {} -c {} -i {} -o {} -r {} -x {} -e {} -inf {}'.format(
                self.production_engine,
                self.system_name + '.prmtop',
                self.system_name + '.inpcrd',
                self.production_prefix + '.{:04d}.rst'.format(window-1),
                self.production_prefix + '.in',
                self.production_prefix + '.{:04d}.out'.format(window),
                self.production_prefix + '.{:04d}.rst'.format(window),
                self.production_prefix + '.{:04d}.nc'.format(window),
                self.production_prefix + '.{:04d}.mden'.format(window),
                self.production_prefix + '.{:04d}.mdinfo'.format(window)
            )
            print('{:<25} {:<10}'.format('Window ' + str(window), self.now()), file=self.log)
            script = open(self.production_prefix + '.{:04d}.sh'.format(window), 'w')
            script.write('#!/usr/bin/env bash\n')
            script.write('source $AMBERHOME/amber.sh\n')
            if self.cuda_visible_devices is not None:
                script.write('export CUDA_VISIBLE_DEVICES={}\n'.format(str(self.cuda_visible_devices)))
            script.write(command)
            script.close()
            # Execute
            sp.call(['bash ./{}.{:04d}.sh'.format(self.production_prefix, window)], shell=True)
            # Check
            print('{:<25} {:<10}'.format('...checking', self.now()), file=self.log)
            try:
                sp.check_output(['grep -q TIMINGS {}.{:04d}.out'.format(self.production_prefix, window)], shell=True)
                print('{:<25} {:<10}'.format('...done', self.now()), file=self.log)
            except sp.CalledProcessError as e:
                print('{:<25} {:<10}'.format('...failed', self.now()), file=self.log)
                sys.exit('Uh oh: looks like production failed.')

    def start_logging(self):
        d_date = dt.datetime.now()
        self.log_time = d_date.strftime("%Y-%m-%d-%I-%M-%S-%p")
        self.log = open(self.log_file + '-' + self.log_time + '.log', 'w', 1)

    def end_logging(self):
        self.log.close()

    def print_parameters(self):
        # Manually putting in the print() statements because otherwise self.__dict__ is unordered
        # and it's hard to find the relevant parameters for each stage of processing. Could alternatively
        # make each stage parameters be in separate dictionaries.
        print('# System parameters', file=self.log)
        dir_path = os.path.dirname(os.path.realpath(__file__))
        host = socket.gethostname()
        print('{:<25} {:<10}'.format(host,dir_path), file=self.log)
        print('{:<25} {:<10}'.format('system', self.system_name), file=self.log)
        print('{:<25} {:<10}'.format('minimize', str(self.minimize)), file=self.log)
        print('{:<25} {:<10}'.format('thermalize', str(self.thermalize)), file=self.log)
        print('{:<25} {:<10}'.format('equilibrate', str(self.equilibrate)), file=self.log)
        print('{:<25} {:<10}'.format('production', str(self.produce)), file=self.log)
        print('{:<25} {:<10}'.format('resume', str(self.resume)), file=self.log)
        print('# File parameters', file=self.log)
        if self.minimize:
            print('{:<25} {:<10}'.format('minimization prefix', self.minimization_prefix), file=self.log)
        if self.thermalize:
            print('{:<25} {:<10}'.format('thermalization prefix', self.thermalization_prefix), file=self.log)
        if self.equilibrate:
            print('{:<25} {:<10}'.format('equilibration prefix', self.equilibration_prefix), file=self.log)
        if self.produce:
            print('{:<25} {:<10}'.format('production prefix', self.production_prefix), file=self.log)
        print('# MD engines', file=self.log)
        if self.minimize:
            print('{:<25} {:<10}'.format('minimization engine', self.minimization_engine), file=self.log)
        if self.thermalize:
            print('{:<25} {:<10}'.format('thermalization engine', self.thermalization_engine), file=self.log)
        if self.equilibrate:
            print('{:<25} {:<10}'.format('equilibration engine', self.equilibration_engine), file=self.log)
        if self.produce:
            print('{:<25} {:<10}'.format('production engine', self.production_engine), file=self.log)
        if self.cuda_visible_devices:
            print('{:<25} {:<10}'.format('CUDA', self.cuda_visible_devices), file=self.log)
            nvidia_file = 'nvidia-smi-{}.log'.format(self.log_time)
            nvidia_log = open(nvidia_file, 'w', 1)
            print('{:<25} {:<10}'.format('...details', nvidia_file), file=self.log)
            sp.call(['nvidia-smi'], stdout=self.log)
            nvidia_log.close()
        if self.minimize:
            print('# Minimization parameters', file=self.log)
            print('{:<25} {:<10}'.format('minimization cycles', self.minimization_cycles), file=self.log)
            print('{:<25} {:<10}'.format('minimization switch', self.minimization_switch), file=self.log)
        if self.thermalize:
            print('# Thermalization parameters', file=self.log)
            print('{:<25} {:<10}'.format('thermalization dt', self.thermalization_dt), file=self.log)
            print('{:<25} {:<10}'.format('thermalization start T', self.thermalization_t0), file=self.log)
            print('{:<25} {:<10}'.format('thermalization end T', self.thermalization_t1), file=self.log)
            print('{:<25} {:<10}'.format('thermalization cut off', self.thermalization_cut), file=self.log)
            print('{:<25} {:<10}'.format('thermalization steps', self.thermalization_steps), file=self.log)
        if self.equilibrate:
            print('# Equilibration parameters', file=self.log)
            print('{:<25} {:<10}'.format('equilibration dt', self.equilibration_dt), file=self.log)
            print('{:<25} {:<10}'.format('equilibration T', self.equilibration_temp), file=self.log)
            print('{:<25} {:<10}'.format('equilibration steps', self.equilibration_steps), file=self.log)
            print('{:<25} {:<10}'.format('equilibration out', self.equilibration_out), file=self.log)
            print('{:<25} {:<10}'.format('equilibration time', self.equilibration_time), file=self.log)
        if self.produce:
            print('# Production parameters', file=self.log)
            print('{:<25} {:<10}'.format('production dt', self.production_dt), file=self.log)
            print('{:<25} {:<10}'.format('production T', self.production_temp), file=self.log)
            print('{:<25} {:<10}'.format('production steps', self.production_steps), file=self.log)
            print('{:<25} {:<10}'.format('production out', self.production_out), file=self.log)
            print('{:<25} {:<10}'.format('production time', self.production_time), file=self.log)
            print('{:<25} {:<10}'.format('production water', str(self.production_water)), file=self.log)

# imin  = 0: Run without minimization
# ntx   = 5: Coordinates will be read
# irest = 1: Restart simulation
# ntpr  = $: Every $ steps, output human-readable energy
# ntwr  = $: Every $ steps, output restart file
# ntwe  = $: Every $ steps, output binary energy and temperature
# ntxo  = 2: NetCDF mdcrd file
# ntwx  = $: Every $ steps, output coordinates to trajectory
# ioutfm= 1: NetCDF trajectory
# iwrap = 1: Coordinates and trajectory are wrapped into box
# ntf   = 2: Bond interactions with H are omitted
# ntc   = 2: SHAKE for hydrogen bonds
# cut   = $: Nonbonded cutoff
# ntt   = 3: Langevin dymamics using gamma_ln (cannonical ensemble)
#       = 0: Constant total energy
#       = 1: Constant temperature (weak coupling)
#       = 2: Anderson coupling (cannoncial ensemble)
#       = 9: Isokinetic Nose-Hoover
# temp0 = $: Reference temperature for ntt > 0
# gamma_ln : Collision frequency for ntt = 3
# ig    =-1: Random number seed
# ntp   = 1: MD with isotropic pressure scaling
#       = 0: No pressure scaling
#       = 2: MD with anisotropic pressure scaling
#       = 3: MD with semi-isotropic pressure scaling
# barostate: Which barostatic
#       = 1: Berendsen
#       = 2: Monte Carlo barostat
# nstlim= $: $ number of MD steps 
# dt    = $: The time step in picoseconds
# nmropt= 1: NMR restraints


a = simulation()
a.simulate()
