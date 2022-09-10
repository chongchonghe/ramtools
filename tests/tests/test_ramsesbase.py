import sys
import ramtools as rt

def main(jobdir):
    
    r = rt.RamsesBase(jobdir)
    outs = r.get_all_outputs()
    times = [r.get_time(i) for i in outs]
    print(f"The list of outputs are: {outs}")
    print(f"The times are {times}")
    print(f"The ds file or the first available output: {r.ds1}")

if __name__ == '__main__':
    print("---------------- Running test_ramsesbase.py ---------------------")
    main(sys.argv[1])
    print('Tests completed successful!')
    print("-----------------------------------------------------------------")
    print()
