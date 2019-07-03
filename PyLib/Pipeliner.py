#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys
import logging
import subprocess
import shlex
import shutil
import time
from inspect import getframeinfo, stack
import threading

logger = logging.getLogger(__name__)


def run_cmd(cmd, ignore_error=False):

    logger.info("Running: " + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error("Error: {}, exit val: {}".format(str(e), e.returncode))

        if ignore_error:
            return e.returncode  # caller decides how to handle the error.
        else:
            raise e

    return 0 # all good.


class Pipeliner(object):

    _checkpoint_dir = None
    _cmds_list = []

    def __init__(self, checkpoint_dir):

        checkpoint_dir = os.path.abspath(checkpoint_dir)

        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
            
        self._checkpoint_dir = checkpoint_dir
    


    def add_commands(self, cmds_list):

        for cmd in cmds_list:
            # check it's a proper Command object
            if not (isinstance(cmd, Command) or isinstance(cmd, ParallelCommandList) ):
                errmsg = "Pipeliner::add_commmands - Error, cmd {} is not a Command or ParallelCommandList object".format(cmd)
                logger.critical(errmsg)
                raise(RuntimeError(errmsg))
            
            self._cmds_list.append(cmd)

    
    def num_cmds(self):
        return len(self._cmds_list)


    def run(self):
        for cmd in self._cmds_list:
            
            checkpoint_dir = self._checkpoint_dir
            cmd.run(checkpoint_dir)

        # since all commands executed successfully, remove them from the current cmds list
        self._cmds_list = list()
    
        return



class Command(object):

    def __init__(self, cmd, checkpoint, ignore_error=False):
        self._cmd = cmd
        self._checkpoint = checkpoint
        self._ignore_error = ignore_error
        self._stacktrace = self._extract_stack(stack())

    def get_cmd(self):
        return self._cmd

    def get_checkpoint(self):
        return self._checkpoint

    def get_ignore_error_setting(self):
        return self._ignore_error
 

    def __repr__(self):
        return self._cmd

    def get_stacktrace(self):
        return self._stacktrace

    def _extract_stack(self, stack_list):
        stacktrace = ""
        for frame_entry in stack_list:
            caller = getframeinfo(frame_entry[0])
            stacktrace += "st: file:{}, lineno:{}\n".format(caller.filename, caller.lineno)

        return stacktrace



    def run(self, checkpoint_dir):

        checkpoint_file = os.path.sep.join([checkpoint_dir, self.get_checkpoint()])
        ret = 0
        if os.path.exists(checkpoint_file):
            logger.info("CMD: " + self.get_cmd() + " already processed. Skipping.")
        else:
            # execute it.  If it succeeds, make the checkpoint file
            start_time = time.time()

            cmdstr = self.get_cmd()
            ret = run_cmd(cmdstr, True)
            if ret:
                # failure occurred.
                errmsg = str("Error, command: [ {} ] failed, stack trace: [ {} ] ".format(cmdstr, self.get_stacktrace()))
                logger.error(errmsg)

                if self.get_ignore_error_setting() is False:
                    raise RuntimeError(errmsg)
            else:
                end_time = time.time()
                runtime_minutes = (end_time - start_time) / 60
                logger.info("Execution Time = {:.2f} minutes. CMD: {}".format(runtime_minutes, cmdstr))
                run_cmd("touch {}".format(checkpoint_file))  # only if succeeds.

        return ret


#############################
## Parallel command execution
#############################


class ParallelCommandThread(threading.Thread):

    def __init__(self, cmdobj, checkpointdir, paraCmdListObj):

        threading.Thread.__init__(self)

        self._cmdobj = cmdobj
        self._checkpointdir = checkpointdir
        self._paraCmdListObj = paraCmdListObj


    def run(self):

        ret = self._cmdobj.run(self._checkpointdir)

        # track job completion and capture success/failure
        self._paraCmdListObj._num_running -= 1
        
        if ret != 0:
            self._paraCmdListObj._num_errors += 1


class ParallelCommandList(object):

    def __init__(self, cmdlist, checkpoint, num_threads, ignore_error=False):

        self._cmdlist = cmdlist
        self._checkpoint = checkpoint
        self._num_threads = num_threads
        self._ignore_error = ignore_error
        self._num_running = 0
        self._num_errors = 0

    def run(self, checkpoint_dir):

        parallel_job_checkpoint_file = self._checkpoint

        full_path_parallel_job_checkpoint_file = os.path.sep.join([checkpoint_dir, parallel_job_checkpoint_file])
        if os.path.exists(full_path_parallel_job_checkpoint_file):
            logger.info("Parallel command series already completed, so skipping. Checkpoint found as: {}".format(full_path_parallel_job_checkpoint_file))
            return
        
        ## run parallel command series, no more than _num_threads simultaneously.
        
        cmd_idx = 0
        while cmd_idx < len(self._cmdlist):

            if self._num_running < self._num_threads:
                
                cmdstr = self._cmdlist[cmd_idx]

                checkpoint_file = "{}.tid-{}".format(parallel_job_checkpoint_file, cmd_idx)

                cmdobj = Command(cmdstr, checkpoint_file, ignore_error=True)
                cmdthread = ParallelCommandThread(cmdobj, checkpoint_dir, self)
                self._num_running += 1
                cmdthread.start() # will auto-decrement _num_threads once it completes via shared memory usage.
                cmd_idx += 1

        ## wait for the remaining ones to finish.
        iter_counter = 0
        while self._num_running > 0:
            time.sleep(1) # avoid unnecessarily running at 100% CPU. Sleeping 1 second is a long enough time in the quantum realm.
            iter_counter += 1
            if iter_counter % 60 == 0:
                sys.stderr.write("\r waiting for {} jobs to complete.    ".format(self._num_running))


        if self._num_errors > 0:
            errmsg = "Error, {} commands failed".format(self._num_errors)
            logger.error(errmsg)
            if not self._ignore_error:
                raise RuntimeError(errmsg)
        else:
            logger.info("All parallel commands succeeded.")
        
        run_cmd("touch {}".format(full_path_parallel_job_checkpoint_file))

        logger.info("done running parallel command series.")
        
        return


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    
    checkpoint_dir = "/tmp/checkpoints_dir." + str(time.time())
    
    pipeliner = Pipeliner(checkpoint_dir)

    pipeliner.add_commands([Command("echo hello!", "hello.ok")])

    pipeliner.add_commands([Command("echo done testing pipeliner", "test.ok")])

    max_x = 10
    cmdlist = [ "echo {} && sleep {}".format(x, max_x + 1 - x) for x in range(1,max_x) ]
    num_threads = 4
    pipeliner.add_commands([ParallelCommandList(cmdlist, "trypara.ok", num_threads)])
    
    
    pipeliner.run()

    shutil.rmtree(checkpoint_dir)

    sys.exit(0)
    
