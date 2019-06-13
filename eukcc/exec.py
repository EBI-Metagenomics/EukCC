import os
import subprocess
from .fileoperations import file

class run():
    '''
    class to handle and run external software
    check if a software exists
    handle input and output
    '''
    def __init__(self, program, inf, outf):
        # check software is in path:
        if run.which(program) is None:
            print("{} is not installed".format(program))
        self.program = program
        self.input = inf
        self.output = outf
        # create output dir
        dc = file.isdir(os.path.dirname(self.output))

    def doIneedTorun(self, force=False):
        if force:
            return(True)
        else:
            return(file.isnewer(self.input, self.output))

    def which(program):
        # taken from
        # https://stackoverflow.com/questions/377017/\
        # test-if-executable-exists-in-python
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

    def run(self, inp, outp, params=[]):
        lst = [self.program]
        lst.extend(params)
        lst.append(inp)
        try:
            with open(outp, "w") as f:
                subprocess.run(lst, stdout=f, check=True)
            return(True)
        except subprocess.CalledProcessError as e:
            print("an error occured while executing {}".format(self.program))
            print(e.output)
            return(False)
        
    def checkOutput(self, file, comment = "#"):
        """
        generic function to see if we have at least one none comment line
        in the output
        """
        with open(file) as f:
            for line in f:
                if not line.startswith(comment):
                    return(True)
        return(False)


# defining classes based on run for executing gmes and hmmscan
class gmes(run):
    def run(self, cores=1):
        lst = [self.program, "-i", self.input,
               "-o", self.output,
               "-n", str(cores)]
        try:
            subprocess.run(lst,  check=True, shell=True)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)


class hmmer(run):
    def run(self, stdoutfile, hmmfiles, cores=1):
        lst = [self.program, "--cpu", str(cores),
               "-o", stdoutfile,
               "--tblout", self.output,
               hmmfiles, self.input]
        #print(" ".join(lst))
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)
    

    
    def clean(self):
        return


def checkDependencies(dep):
    """
    simple function to check all dependencies
    """
    allgood = True
    for name, program in dep.items():
        if run.which(program) is None:
            print("Dependency: Could not find {} ({})".format(program, name))
            allgood = False

    return(allgood)
