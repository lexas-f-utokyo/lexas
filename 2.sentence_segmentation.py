import scispacy
import en_core_sci_sm
en = en_core_sci_sm.load()


f = open("data/result_sections_segmentation.txt", "w")
f2 = open("data/result_sections.txt", "r")


def segmentation(sent):
    doc = en(sent)
    ret = "#####".join([str(a) for a in doc.sents])
    return ret


if __name__ == "__main__":
    for line in f2:
        ls = line.strip().split("\t")
        if len(ls) == 3:
            sent = ls[2]
            if len(sent) > 1000000:
                sent = sent[:1000000]
            seg = segmentation(sent)
            f.write(ls[0] + "\t" + ls[1] + "\t" + seg + "\n")
f.close()
f2.close()
