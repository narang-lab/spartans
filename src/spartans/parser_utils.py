def _spartans_logo():
    '''
    ascii logo from https://patorjk.com/software/taag/
    '''
    return r"""
______             ____ _____      _   _______
\  ___)           |  _ (_   _)    | \ | \  ___)
 \ \  ______ __  _| |_) )| | __  _|  \| |\ \
  > >(  __  )  \/ /  __/ | |/  \/ /     | > >
 / /__| || ( ()  <| |    | ( ()  <| |\  |/ /__
/_____)_||_|\__/\_\_|    |_|\__/\_\_| \_/_____)

"""

def print_spartans_logo(log):
    log.info(_spartans_logo())
