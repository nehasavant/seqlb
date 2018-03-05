#!/usr/bin/env python

import seqlib

##Maybe add a parse command line? 
def parse_command_line():
    " parses args for the helloworld.py funtions"

    # init parser and add arguments
    parser = argparse.ArgumentParser()

    # add short args
    parser.add_argument(
        "-n", "--name",
        help="optional name to say hello to",
        type=str)

def main():
    " run seqlib library"
    args = parse_command_line()
    seqlib.seqlib(
        name=args.name,
        howlong=args.howlong,
        countdown=args.countdown)
	    