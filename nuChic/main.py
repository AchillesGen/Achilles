def main():
    ''' 
    The Main prgram code
    Authors:

    '''
    logging.info('Starting nuChic')

    logging.info('nuChic ended succesfully')

if __name__ == '__main__':
    # initialize the logger to keep track of what happens
    import logging.config
    logger_dict = {
            'version': 1,
            'disable_existing_loggers': False,
            'loggers': {
                '': {
                    'level': 'DEBUG',
                    },
                },
            }

    logging.config.dictConfig(logger_dict)

    main()
