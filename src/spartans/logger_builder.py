import logging

def setup_logger(logger_name, log_file, level=logging.INFO):
        '''
        Factory method for building logging objects. Code from:
        http://stackoverflow.com/questions/17035077/python-logging-to-multiple-log-files-from-different-classes
        '''

        l = logging.getLogger(logger_name)
        formatter = logging.Formatter('%(asctime)s : %(message)s')
        fileHandler = logging.FileHandler(log_file, mode='w')
        fileHandler.setFormatter(formatter)
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)

        l.setLevel(level)
        l.addHandler(fileHandler)
        l.addHandler(streamHandler)

        # return l
        # Note: we don't return logger object to reduce risk of this warning from the documentation:
        # Loggers should NEVER be instantiated directly, but always through the module-level function logging.getLogger(name).

