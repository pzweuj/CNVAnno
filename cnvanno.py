# coding=utf-8
# pzw
# 20241213

import sys
import argparse
from src.utils import cnv_anno_process
from src.database_prepare import prepare_all

# 运行
def main():
    parser = argparse.ArgumentParser(
        description="Usage：python3 cnvanno.py anno -i <input bed> -o <output txt> -g <hg19|hg38>",
        prog="cnvanno.py",
        usage="python3 cnvanno.py -h",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.1 20241213")
    parser.add_argument('option', type=str,
                        choices=["anno", "database"],
                        help='\nOptions：\nanno      annotation mode\ndatabase     database prepare mode\n')
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args(sys.argv[1:2])
    sys.argv = sys.argv[:1] + sys.argv[2:]
    option_set = args.option
    
    if option_set == "anno":
        parser.add_argument("-i", "--input", type=str, help="input bed，format: chr1 100 200 DUP")
        parser.add_argument("-o", "--output", type=str, help="output txt")
        parser.add_argument("-g", "--genome", type=str, help="genome version，default is hg19，[hg19, hg38]", default="hg19")
    elif option_set == "database":
        parser.add_argument("-g", "--genome", type=str, help="genome version，[hg19, hg38]")
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    
    args = parser.parse_args()
    if option_set == "anno":
        cnv_anno_process(args.input, args.output, args.genome)
    elif option_set == "database":
        prepare_all(args.genome)

if __name__ == "__main__":
    main()

