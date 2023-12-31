import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data")
    args = parser.parse_args()

    data = pd.read_csv(args.data, header=None, names=["score", "x", "y"])

    x = data["x"]
    x_max = np.max(x) // 10 + 1

    y = data["y"]
    y_max = np.max(y) // 10 + 1

    print(x_max, y_max)

    colormap = np.zeros((x_max, y_max), dtype=np.float32)

    for s, x, y in zip(data["score"], data["x"], data["y"]):
        colormap[int(x // 10), int(y // 10)] = s

    print(colormap)

    im1 = plt.imshow(colormap.T)
    plt.colorbar(im1)
    plt.show()


if __name__ == "__main__":
    main()
