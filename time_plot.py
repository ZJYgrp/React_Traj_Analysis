import numpy as np
from matplotlib import pyplot


def plot_gaps(times):
    x_formation_gap = list(map(lambda result:
                               abs(result[2] - result[1]) if result[0] == 'X' else None, times))
    y_formation_gap = list(map(lambda result:
                               abs(result[2] - result[1]) if result[0] == 'Y' else None, times))

    bins = np.arange(min(min(x_formation_gap), min(y_formation_gap)),
                     max(max(x_formation_gap), max(y_formation_gap)) + 5, 5)
    pyplot.hist(x_formation_gap, bins, alpha=0.5, label='X')
    pyplot.hist(y_formation_gap, bins, alpha=0.5, label='Y')
    pyplot.legend(loc='upper right')
    pyplot.axvline(x=np.median(x_formation_gap), color='b', linestyle='dashed', linewidth=2)
    pyplot.axvline(x=np.median(y_formation_gap), color='orange', linestyle='dashed', linewidth=2)
    # pyplot.title('Bond Formation Time Gaps')
    # pyplot.xlabel('Time Gap (fs)')
    # pyplot.ylabel('Frequency')
    # pyplot.show()
    pyplot.savefig('./trajTS/bond_time_gaps.png', dpi=100)


def plot_bonds(times):
    x_formation_time = list(map(lambda result: result[1] if result[0] == 'X' else None, times))
    y_formation_time = list(map(lambda result: result[1] if result[0] == 'Y' else None, times))

    bins = np.arange(min(x_formation_time), max(y_formation_time) + 5, 5)

    pyplot.hist(x_formation_time, bins, alpha=0.5, label='X')
    pyplot.hist(y_formation_time, bins, alpha=0.5, label='Y')
    pyplot.legend(loc='upper right')
    pyplot.axvline(x=np.median(x_formation_time), color='b', linestyle='dashed', linewidth=2)
    pyplot.axvline(x=np.median(y_formation_time), color='orange', linestyle='dashed', linewidth=2)
    # pyplot.title('Bond Formation Times')
    # pyplot.xlabel('Formation Time (fs)')
    # pyplot.ylabel('Frequency')
    # pyplot.show()
    pyplot.savefig('./trajTS/bond_times.png', dpi=100)


# main func
def main():
    time_file = open('./trajTS/trajTimes.txt')
    times = list(map(lambda x:
                     list(map(lambda y:
                              float(y) if y.replace('.', '', 1).isdigit() else y,
                              x.split(','))),
                     time_file.readlines()))
    print(times)
    time_file.close()

    if len(times[0]) == 3:
        plot_gaps(times)
    elif len(times[0]) == 2:
        plot_bonds(times)
    time_file.close()


if __name__ == '__main__':
    main()
