import logging
logger = logging.getLogger('main') # use logger instead of printing


def main_func(params): #with a name that is related to the file's name
    pass # do something

def auxiliary_func1(params):
    pass # do something

def auxiliary_func2(params):
    pass # do something

#etc...


if __name__ == '__main__':
    from sys import argv
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    if len(argv) < 4: # change the number of arguments according to your script
        logger.error('Usage: python ' + argv[0] + ' <arg1> <arg2> <arg3>...')
        exit()
    else:
        main_func(*argv[1:])