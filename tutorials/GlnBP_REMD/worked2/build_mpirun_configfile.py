#!/usr/bin/env python
import os
import re
import argparse
import textwrap
import abc
import subprocess as sp
import collections
import warnings

ABC = abc.ABCMeta('ABC', (object,), {})  # compatible with Python 2 *and* 3

warnings.filterwarnings("once")

correct_major_mpich = 3


def check_mpich():
    """
    Try to check for MPICH3 installation over MPICH2.
    Because MPICH is not required to run this module, we do not enforce its installation.
    However, this moodule only produces MPICH3 configfiles, so we warn the users if there is a problem
    """

    auto_detect_warning = "Unable to auto-detect MPICH Version.\n" \
                          "Your host and configfiles creation will still be attempted, but note that " \
                          "build_mpirun_configfile only builds MPICH{correct_major} compatible files."
    wrong_major_warning = "Detected MPICH version {wrong_major}! \n" \
                          "Your host and configfiles creation will still be attempted, but you may have problems as " \
                          "build_mpirun_configfile only builds MPICH{correct_major} compatible files."

    try:
        # Send the error output (if any) to null. Because this check is 100% optional, we dont want un-needed errors
        # showing up in people's logs
        with open(os.devnull, 'w') as devnull:  # Python 2 and 3 compatible
            shell_output = sp.check_output("mpichversion", shell=True, stderr=devnull)
    except sp.CalledProcessError:
        try:
            # Try the mpich2 fallback
            with open(os.devnull, 'w') as devnull:
                discarded = sp.check_output("mpich2version", shell=True, stderr=devnull)
            print(wrong_major_warning.format(wrong_major=2, correct_major=correct_major_mpich))
            return
        except sp.CalledProcessError:
            print(auto_detect_warning.format(correct_major=correct_major_mpich))
            return
    txt_output = bytestring_to_string(shell_output)
    regex_match = re.search("Version[^\d]*(\d)", txt_output)
    if regex_match is None:
        print(auto_detect_warning.format(correct_major=correct_major_mpich))
        return
    major_version = int(regex_match.group(1))
    if major_version != correct_major_mpich:
        print(wrong_major_warning.format(wrong_major=major_version, correct_major=correct_major_mpich))
    return


def bytestring_to_string(input_string):
    """Convert byte strings to python string, e.g. from subprocess"""
    return input_string.decode("utf-8").strip()


class HydraConfigFileCreator(ABC):
    def __init__(self, mpiversion):
        """
        Parent class handling all the Hydra config file creation
        Sublcasses will be the file system specific handlers
        """
        self.cuda_visible_devices = collections.OrderedDict()
        self.hydra_delimiter = self._set_hydra_delimiter(mpiversion)
        return

    def build_hydra_configfile_hosts(self, host_list=None):
        """Return the hosts line of the configfile, given a Python list of host names"""
        if host_list is None:
            host_list = self.extract_hostlist
        output_string = "-hosts "
        output_string += ":1,".join(host_list) + ":1"
        return output_string

    def _set_hydra_delimiter(self, mpiversion):
        if mpiversion == '3':
            delimiter = '\n'
        else:
            delimiter = ':\n'
        return delimiter

    def write_configfile(self, command_list, configfile_output_filepath='configfile',
                         hostfile_output_filepath='hostfile'):
        # Add Hosts line
        base_conf_string = "-np 1 -env CUDA_VISIBLE_DEVICES {cvd} {cmd} {delimiter}"
        base_host_string = "{host}\n"
        proc_list = []
        host_list = []
        command_str = ' '.join(command_list)
        for host, cvd in self.extract_host_cuda_visible_devs():
            conf_string = base_conf_string.format(cvd=cvd, cmd=command_str, delimiter=self.hydra_delimiter)
            host_string = base_host_string.format(host=host)
            proc_list.append(conf_string)
            host_list.append(host_string)
        # Strip last colon from the gpu_proc list
        proc_list[-1] = proc_list[-1][:-len(self.hydra_delimiter)]
        config_output_string = "".join(proc_list) + "\n\n"  # Extra \n since delimiter was removed
        hostfile_output_string = "".join(host_list) + "\n"
        with open(configfile_output_filepath, 'w') as conf_out_file:
            conf_out_file.write(config_output_string)
        with open(hostfile_output_filepath, 'w') as host_out_file:
            host_out_file.write(hostfile_output_string)

    @abc.abstractmethod
    def extract_hostlist(self):
        """
        Should return a Python list of strings of hosts from the hostfile
        Each host should be replicated for each process its in
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def extract_host_cuda_visible_devs(self):
        """
        Returns a list of [host, CUDA_VISIBLE_DEVICE] lists where the last entry is a single integer
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def hostfile(self):
        """Should return the string that is the host list on file"""
        raise NotImplementedError()


class PBSHydraConfigCreator(HydraConfigFileCreator):

    def __init__(self, *args):
        with open(os.environ.get('PBS_GPUFILE'), 'r') as gpu_file:
            self._hostfile = gpu_file.read()
        # $PBS_GPUFILE format is like 'gpu-1-13-gpu0', where '%(hostname)s-gpu%(gpuid)s'
        lines = self.hostfile.split('\n')
        hosts = []
        gpus = []
        for line in lines:
            if line == "":
                continue
            line.strip()
            result = re.match('^(\S+)-gpu(\d)$', line)
            hosts.append(result.group(1))
            gpus.append(result.group(2))
        self._host_list = hosts
        self._gpu_list = gpus
        super(PBSHydraConfigCreator, self).__init__(*args)

    def extract_hostlist(self):
        return self._host_list

    def extract_host_cuda_visible_devs(self):
        # Sift through the raw files
        cvd_list = []
        for host, gpu in zip(self._host_list, self._gpu_list):
            cvd_list.append([host, gpu])
        return cvd_list

    @property
    def hostfile(self):
        return self._hostfile


class LSFHydraConfigCreator(HydraConfigFileCreator):

    def extract_host_cuda_visible_devs(self):
        host_list = self.hostfile.split()
        host_set = set(host_list)
        host_count = collections.Counter(host_list)
        mps_set = collections.OrderedDict()
        for host in host_set:
            # Get CUDA_VISIBLE_DEVICES from each host
            # Do not reverse the order of the quotes here!
            # The remote command MUST be in single quotes or the host will convert the string to its own local variable
            # then tell the remote to echo a literal string instead of fetching its own version of the variable
            cvdh = sp.check_output("blaunch -z {} 'echo $CUDA_VISIBLE_DEVICES'".format(host), shell=True)
            # Cast into string
            cvdh = bytestring_to_string(cvdh)
            self.cuda_visible_devices[host] = cvdh
            # Check on MPS
            mps = sp.check_output("blaunch -z {} 'echo $LSB_START_JOB_MPS'".format(host), shell=True)
            mps_set[host] = bytestring_to_string(mps)
        # Cast CUDA_VISIBLE_DEVICES by host
        cvd_list = []
        for host in host_set:
            procs_on_host = host_count[host]
            if mps_set[host] == "Y":
                looped_cvd = [i for i in range(procs_on_host)]  # Force cast to list
            else:
                looped_cvd = [int(i) for i in self.cuda_visible_devices[host].split(',')]
            for n in range(procs_on_host):
                cvd_list.append([host, looped_cvd[n]])
        return cvd_list

    def extract_hostlist(self):
        return self.hostfile.split()

    @property
    def hostfile(self):
        return os.environ.get('LSB_HOSTS')

class SLURMHydraConfigCreator(HydraConfigFileCreator):
  
    def __init__(self, *args):
        super(SLURMHydraConfigCreator, self).__init__(*args)
        # Get out the nodelist and CVD for this proces
        # This requires doing srun over python so that the local environment variables are not processed first
        # This long string avoids creating an extra file that `srun` can process
        # import os, execute 2 os.environ commands, stitch them together with a space.
        # the "+" are inside the print command to make the 3 strings come together.
        hosts_cvd = sp.check_output("srun python -c 'import os; print(os.environ.get(\"SLURMD_NODENAME\") + "
                                    "\" \" + "
                                    "os.environ.get(\"CUDA_VISIBLE_DEVICES\"))'", shell=True)
        hosts_cvd = bytestring_to_string(hosts_cvd).split('\n')
        host_list = []
        for host_cvd_pair in hosts_cvd:
            host, cvd = host_cvd_pair.split(' ')
            host_list.append(host)
            self.cuda_visible_devices[host] = cvd
        self._host_list = tuple(host_list)

    def extract_host_cuda_visible_devs(self):
        host_list = self.extract_hostlist()
        # Ensure we get an ordered "set", normal set is un-ordered
        # We don't care that its a dict, we only need the keys
        # This fixes a rare bug appearing to involve SLURM and Hydra MPI's --control-port when the first node
        # in the config file is not the control-port
        host_set = collections.OrderedDict.fromkeys(sorted(host_list))
        cvd_list = []
        for host in host_set:
            cvds = self.cuda_visible_devices[host]
            for cvd in cvds.split(','):
                cvd_list.append([host, cvd])
        return cvd_list

    def extract_hostlist(self):
        return self._host_list

    @property
    def hostfile(self):
        return os.environ.get('SLURM_JOB_NODELIST')


def main():
    args, exec_args = parse_args()
    manager = figure_out_manager(args.mpiversion)
    manager.write_configfile(exec_args, configfile_output_filepath=args.configfilepath,
                             hostfile_output_filepath=args.hostfilepath)


def figure_out_manager(mpiversion):
    """Figure out what manager we are using by trial and error"""
    # PBS
    pbs_gpu_file = os.environ.get('PBS_GPUFILE')
    if pbs_gpu_file is not None:
        return PBSHydraConfigCreator(mpiversion)
    lsf_host_file = os.environ.get('LSB_HOSTS')
    if lsf_host_file is not None:
        return LSFHydraConfigCreator(mpiversion)
    slurm_host_file = os.environ.get('SLURM_JOB_NODELIST')
    if slurm_host_file is not None:
        return SLURMHydraConfigCreator(mpiversion)
    raise RuntimeError("Cannot determine job scheduler!\n"
                       "Please ensure one of the following environment variables is set for your job:\n"
                       "    PBS: \"PBS_GPUFILE\"\n"
                       "    LSF: \"LSB_HOSTS\"\n"
                       "  SLURM: \"SLURM_JOB_NODELIST\"")


def parse_args():
    help_text = textwrap.dedent(r"""
        Construct a configfile for MPICH3 mpirun from one of the following:
            - Torque/Moab $PBS_GPUFILE
            - LSF $LSB_HOSTS

        Put the command to be executed after the command for this script, e.g.
        "build_mpirun_configfile python yourscript.py -yourarg x".
        Run in a batch job script or interactive session as follows:
            python build_mpirun_configfile.py python yourscript.py -yourarg x
            mpirun -f hostfile -configfile configfile"""
    )
    argparser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)
    argparser.add_argument(
        '--mpiversion', type=str, choices=[correct_major_mpich], default=correct_major_mpich,
        help="MPICH2 and MPICH3 have different line wrap specifications, "
             "This setting allows you to control which version of wrapping to use. "
             "MPICH2 IS NO LONGER SUPPORTED"
             "(default: 3)"
    )
    argparser.add_argument(
        '--configfilepath', type=str, default='configfile',
        help='optionally specify a filepath for the configfile (default: "configfile")'
    )
    argparser.add_argument(
        '--hostfilepath', type=str, default='hostfile',
        help='optionally specify a filepath for the hostfile (default: "hostfile")'
    )
    argparser.add_argument(
        '--nocheckmpich', action='store_false',
        help='optionally dont check for mpich version'
    )
    args, unknown_args = argparser.parse_known_args()
    if args.nocheckmpich:  # Setting this arg makes it false
        check_mpich()
    exec_args = unknown_args
    return args, exec_args

if __name__ == "__main__":
    main()
