from nudel import Nuclide

if __name__ == '__main__':

    nuc = Nuclide(48,23)
    print(nuc)
    print(nuc.mass)
    print(nuc.protons)
    #print(nuc.adopted_levels.records)
    print(nuc.adopted_levels.levels)
    print(list(nuc.get_isomers()))
    print(list(nuc.get_daughters()))