from toolshed import reader, nopen
import sys

def stream_file(f, pad_info=None):
    if pad_info is None:
        pad_info = dict(upstream=0, downstream=0)

    if pad_info['upstream'] == pad_info['downstream'] == 0:
        for line in nopen(f):
            yield line.rstrip("\r\n")
        raise StopIteration

    for toks in reader(f, header=None):
        toks[1] = int(toks[1])
        toks[2] = int(toks[2])

        if pad_info['upstream'] != 0:
            tok_strand = '-' if len(toks) > 3 and (toks[3] == '-' or toks[-1] == '-') else '+'

            if tok_strand == '+':
                toks[1] += pad_info['upstream']
            else:
                toks[2] -= pad_info['upstream']

        if pad_info['downstream'] != 0:
            tok_strand = '-' if len(toks) > 3 and (toks[3] == '-' or toks[-1] == '-') else '+'

            if tok_strand == '+':
                toks[2] += pad_info['downstream']
            else:
                toks[1] -= pad_info['downstream']

        toks[1], toks[2] = str(max(0, toks[1])), str(max(0, toks[2]))
        yield "\t".join(toks)


def parse_file_pad(fname):
    """
    >>> parse_file_pad("asdf.bed:100")
    {'downstream': 100, 'file': 'asdf.bed', 'upstream': -100}

    >>> parse_file_pad("asdf.bed:100:100")
    {'downstream': 100, 'file': 'asdf.bed', 'upstream': -100}

    >>> parse_file_pad("asdf.bed:-100:100")
    {'downstream': 100, 'file': 'asdf.bed', 'upstream': 100}

    >>> parse_file_pad("asdf.bed:-100:0")
    {'downstream': 0, 'file': 'asdf.bed', 'upstream': 100}

    >>> parse_file_pad("asdf.bed:0:-100")
    {'downstream': -100, 'file': 'asdf.bed', 'upstream': 0}

    """

    def parse_dist(dist):
        strand = None
        if dist[0] in "+-":
            strand = dist[0]
            dist = dist[1:]
        return int(dist), strand

    toks = fname.split(":")
    if len(toks) > 3:
        1/0

    f, upstream, downstream = toks[0], 0, 0

    if len(toks) > 1:
        upstream, upstrand = parse_dist(toks[1])
        downstream = -upstream if upstrand == "-" else upstream
        upstream = upstream if upstrand == "-" else -upstream

    if len(toks) > 2:
        downstream, downstrand = parse_dist(toks[2])
        downstream = -downstream if downstrand == "-" else downstream

    return dict(file=f,
            upstream=upstream,
            downstream=downstream)

if __name__ == "__main__":
    import doctest
    doctest.testmod()


    info = parse_file_pad(sys.argv[1])
    print info["file"], info

    for s in stream_file(info['file'], info):
        print s
