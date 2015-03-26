import time

class Timer(object):
    """ Timer object that can be used
    in the with context.

    with Timer('my timer'):
        do_expensive_operation()

    from: http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
    """

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)
