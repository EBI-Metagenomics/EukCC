import os
import logging
import shutil
from eukcc.base import which
from eukcc.file import file
from eukcc.fasta import clean_metaeuk_fasta
import subprocess


class run:
    """
    class to handle and run external software
    """

    def __init__(
        self, program, workdir, infiles, outfiles, cores=1, touch=False, **kwargs
    ):
        # check software is in path:
        if which(program) is None:
            raise EnvironmentError("Could not find executable: {}".format(program))
        self.success = None

        if not isinstance(infiles, list):
            infiles = [infiles]
        if not isinstance(outfiles, list):
            outfiles = [outfiles]

        # create workdir
        file.isdir(workdir)

        # check for all infiles
        for infile in infiles:
            if infile is None or (
                not os.path.exists(infile)
                and not os.path.exists("{}.h3f".format(infile))
            ):
                raise OSError("File not found: {}".format(infile))
        # make sure all infiles are abspath
        infiles = [os.path.abspath(x) for x in infiles]

        # check if running is required
        if self.run_needed(workdir, infiles, outfiles):
            logging.debug("Running {}".format(program))
            self.run(program, workdir, infiles, outfiles, cores, kwargs)

    def touch(self, files=None):
        """
        simple function, to only touch output files
        """
        logging.debug("Only touching now")
        if files is None:
            path = self.output
            with open(path, "a"):
                os.utime(path, None)
        else:
            for path in files:
                with open(path, "a"):
                    os.utime(path, None)

    def run_needed(self, workdir, infiles, outfiles, force=False):
        logging.debug("Testing if I need to run this step")
        if force:
            return True
        else:
            for inf in infiles:
                for p in outfiles:
                    ouf = os.path.join(workdir, p)
                    if not os.path.exists(inf) and os.path.exists("{}.h3f".format(inf)):
                        inf = "{}.h3f".format(inf)
                    _need_run = file.isnewer(inf, ouf)
                    if _need_run:
                        if _need_run > 1:
                            logging.debug(
                                "Need to run because {} does not exists yet".format(ouf)
                            )
                        else:
                            logging.debug(
                                "Need to run because {} is older than {}".format(
                                    ouf, inf
                                )
                            )
                        return True
            return False

    def run(self, program, workdir, infiles, outfiles, cores, **kwargs):
        raise NotImplementedError("needs to be implemented downstream")

    def cleanup(self, workdir, folders=None):
        """Remove temporary folders

        Given a list of folder names we will remove one after
        another, so that temporary directories are gone.
        """

        if folders is not None:
            if type(folders) == str:
                folders = [folders]
            for p in folders:
                path = os.path.join(workdir, p)
                if os.path.exists(path):
                    shutil.rmtree(path)
                    logging.debug("Removed temp folder: %s", path)
                else:
                    logging.debug("Cant remove folder, as it does not exists: %s", path)


class hmmsearch(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "hmmsearch.stdout")
        lst = [
            program,
            "-o",
            "/dev/null",
            "--noali",
            "--cpu",
            str(cores),
            "--tblout",
            outfiles[0],
        ]
        if "cut_ga" in kwargs.keys() and kwargs["cut_ga"] is False:
            lst.append("")
        else:
            lst.append("--cut_ga")
        # add last cmd line args
        lst.extend(
            [
                infiles[0],
                infiles[1],
            ]
        )
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
                self.success = True
        except subprocess.CalledProcessError:
            logging.info("Hmmsearch failed, check logfile {}".format(logfile))
            self.success = False
            exit(1)


class hmmalign(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "hmmalign.stdout")
        lst = [
            program,
            "--amino",
            "--outformat",
            "afa",
            "-o",
            outfiles[0],
            "--mapali",
            infiles[0],
            infiles[1],
            infiles[2],
        ]
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            self.success = True
        except subprocess.CalledProcessError:
            self.success = False
            logging.info("Hmmalign failed, check logfile {}".format(logfile))
            exit(1)


class hmmfetch(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "hmmfetch.stdout")
        lst = [program, "-o", outfiles[0], "-f", infiles[0], infiles[1]]
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            self.success = True
        except subprocess.CalledProcessError:
            logging.info("Hmmfetch failed, check logfile {}".format(logfile))
            self.success = False
            exit(1)


class hmmpress(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "hmmpress.stdout")
        lst = [program, "-f", infiles[0]]
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            self.success = True
        except subprocess.CalledProcessError:
            logging.info("Hmmpress failed, check logfile {}".format(logfile))
            self.success = False
            exit(1)


class epa_split(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "epa-ng_split.stout")
        lst = [program, "--split", infiles[0], infiles[1]]
        # make sure to remove old files first
        for f in outfiles:
            file.remove(os.path.join(workdir, f))
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            self.success = True
        except subprocess.CalledProcessError:
            self.success = False
            logging.info("epa-ng failed, check logfile {}".format(logfile))


class epa_ng(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "epa-ng.stout")
        lst = [
            program,
            "--threads",
            str(cores),
            "--ref-msa",
            infiles[0],
            "--tree",
            infiles[1],
            "--query",
            infiles[2],
            "--model",
            kwargs["model"],
            "--preserve-rooting",
            "on",
            "--redo",
        ]
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            self.success = True
        except subprocess.CalledProcessError:
            self.success = False
            logging.info("epa-ng failed, check logfile {}".format(logfile))


class guppy(run):
    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "guppy.stout")
        lst = [program, "tog", "-o", outfiles[0], infiles[0]]
        try:
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
        except subprocess.CalledProcessError:
            logging.info("guppy tog failed, check logfile {}".format(logfile))
            exit(1)


class metaeuk(run):
    """
    class to run and predict proteins using metaeuk
    Will run  easy-predict
    """

    def run(self, program, workdir, infiles, outfiles, cores=1, kwargs=None):
        logfile = os.path.join(workdir, "metaeuk.stout")
        lst = [
            program,
            "easy-predict",
            "--protein",
            "1",
            "--translation-table",
            "1",
            "--threads",
            str(cores),
            infiles[0],
            infiles[1],
            kwargs["prefix"],
            "metaeuk_tmp",
        ]
        try:
            # remove prediction if it was done before update to the database or MAG
            if file.exists(os.path.join(workdir, outfiles[0])):
                os.remove(os.path.join(workdir, outfiles[0]))
            with open(logfile, "a") as fout:
                subprocess.run(
                    " ".join(lst),
                    cwd=workdir,
                    check=True,
                    shell=True,
                    stdout=fout,
                    stderr=fout,
                )
            # clean the output for epa-ng
            clean_metaeuk_fasta(
                os.path.join(workdir, outfiles[0]), os.path.join(workdir, outfiles[1])
            )
        except subprocess.CalledProcessError:
            logging.info("metaeuk failed, check logfile {}".format(logfile))
            exit(1)
