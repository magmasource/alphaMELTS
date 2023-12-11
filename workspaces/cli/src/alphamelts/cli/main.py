from rich import print

from alphamelts.lib import Foo


def main() -> None:
    print(Foo(name="cli"))


if __name__ == "__main__":
    main()
