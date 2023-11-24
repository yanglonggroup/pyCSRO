import os
from argparse import ArgumentParser, RawTextHelpFormatter

def main():
    _BANNER = """
    This is the Yang Group python packaging format.
    
    This helps people set up their git repo correctly for releasing.
    """
    parser = ArgumentParser(prog = 'pkgname',
                            description = _BANNER, formatter_class = RawTextHelpFormatter)
    # optional args
    parser.add_argument("--opt", dest="opt", type = str,
                        help = "option")

if __name__ == "__main__":
    main()