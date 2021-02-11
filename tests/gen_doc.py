import re

if __name__ == "__main__":
    # parse units.h
    category_map = {}
    literal_map = {}
    categories = []
    units = []
    relationships = []
    with open("../units.h", "r") as fp:
        for line in fp:
            # unit
            match = re.match(R'^using\s*([\w]*)\s*=\s*Unit<([\w]*)', line)
            if match:
                name = match[1]
                category = match[2]
                category_map[name] = category
                units.append(name)
                if(category not in categories):
                    categories.append(category)
            # GEN_LITERAL
            match = re.match(R'^GEN_LITERAL\s*\(\s*(\w*)\s*,\s*(\w*)\s*\)', line)
            if match:
                literal_map[match[2]] = match[1]
            # GEN_INVERSE
            match = re.match(R'^GEN_INVERSE\s*\(\s*(\w*)\s*,\s*(\w*)\s*\)', line)
            if match:
                relationships.append("%s = inverse(%s)" % (category_map[match[1]], category_map[match[2]]))
            # GEN_MULT_DIV
            match = re.match(R'^GEN_MULT_DIV\s*\(\s*(\w*)\s*,\s*(\w*)\s*,\s*(\w*)\s*\)', line)
            if match:
                relationships.append("%s = %s * %s" % (category_map[match[1]], category_map[match[2]], category_map[match[3]]))
            match = re.match(R'^GEN_MULT_DIV_SQ\s*\(\s*(\w*)\s*,\s*(\w*)\s*\)', line)
            if match:
                relationships.append("%s = %s * %s" % (category_map[match[1]], category_map[match[2]], category_map[match[2]]))

    print("num units: ", len(units))
    print("num relationships: ", len(relationships))

    # output units.md
    with open("../units.md","w") as fp:
        fp.write("## Units (and literal)\n\n")
        # write categories and units
        for category in categories:
            fp.write("### " + category + "\n")
            for unit in units:
                if category_map[unit] == category:
                    fp.write("* " + unit)
                    if unit in literal_map:
                        fp.write(" (%s)" % literal_map[unit])
                    fp.write("\n")
            fp.write("\n")
        fp.write("## Relationships\n### Note: multiplies have corresponding divides\n")
        for relationship in relationships:
            fp.write("* " + relationship + "\n")



