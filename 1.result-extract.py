import xml.etree.ElementTree as ET


def parse(pmcid):
    try:
        tree = ET.parse("articles/" + pmcid + ".nxml")

    except FileNotFoundError:
        return "Year", "NotFound"
    except BaseException:
        return "Year", "Error"

    try:
        root = tree.getroot()
        art = root[0].find("article-meta")
        pubd = art.find("pub-date")
        year = pubd.find("year")
        year = year.text
        sent = []
    except BaseException:
        return "Year", "YearError"

    if len(root) < 2:
        return year, "FormatError"

    for sec in root[1]:
        # root[0]:front root[1]:body
        # sec1 : Introduction
        title = sec.find('title')
        if "result" in title.text.lower():
            for child in sec:
                if child.tag == "p":
                    for fig in child.findall('fig'):
                        child.remove(fig)
                    sent.append(" ".join(child.itertext()).replace(
                        "\n", "").replace("\t", " "))
                else:
                    for chi_chi in child:
                        if chi_chi.tag == "p":
                            for fig in chi_chi.findall('fig'):
                                chi_chi.remove(fig)
                            sent.append(" ".join(chi_chi.itertext()).replace(
                                "\n", "").replace("\t", " "))
    return year, " ".join(sent)


f = open("data/result_sections.txt", "w")
with open("pmcid_list.txt") as f2:
    for line in f2:
        pmcid = line.strip()
        year, sent = parse(pmcid)
        f.write("\t".join([year, pmcid, sent]) + "\n")
