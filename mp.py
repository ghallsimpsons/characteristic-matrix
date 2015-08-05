from multiprocessing import cpu_count
import pp
import time

def mpexec(kernel, param, values, **kwargs):
    """
    Execute the kernel on the local machine, splitting values of
    "param" across multiple processes.
    """
    if param in kwargs.keys():
        raise RuntimeError(
            "Sweep parameter '{}' can not be passed as a kwarg.".format(
            param))
    nprocs = cpu_count()-1
    job_server = pp.Server()
    job_server.set_ncpus(nprocs)
    print "Pooling {} procs.".format(job_server.get_ncpus())
    start_time = time.time()

    # Create jobs from parameter space
    jobs = []
    for v in values:
        ext_kwargs = kwargs.copy()
        ext_kwargs.update({param: v})
        def job(kernel, kwargs):
            return kernel(**kwargs)
        jobs.append(job_server.submit(job, (kernel, ext_kwargs)))
    results = [result() for result in jobs]
    end_time = time.time()
    print "elapsed: {}".format((end_time-start_time))
    return results
