import gc


def dump_garbage():
    print "\nGARBAGE"
    gc.collect()

    print "\nGARBAGE OBJECTS:"
    for x in gc.garbage:
        if type(x) is tuple:
            print str(x)
            print "Cell value is = " + str(x[0].cell_contents)
            # print "Function value is =" + str(x[1])
        # s = str(x)
        # if len(s) > 80: s = s[:80]
        # print type(x), "\n ", s


if __name__ == "__main__":
    import gc
    import objective_geoconprecheck as obj
    gc.enable()
    gc.set_debug(gc.DEBUG_LEAK)

    obj.main()

    dump_garbage()