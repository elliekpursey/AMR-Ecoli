#!/usr/bin/env python3

import pandas as pd
import numpy as np

bs_raw = pd.read_csv('../results/biosample_attributes_EC.csv', index_col=0) 

# make lists of defining words to include for each category:
    
# clinical/human/hospital associated

human = ['(Homo)', 'Homo', 'sapiens', 'Human', '(Human)']

human_regex =  ('|'.join(human))

# other wild animal        
ani = ['(Panda)', 'Panda', '(panda)', '(Callithrix)', 'Deer', '(deer)', 'deer', '(Guinea)', '(Daphnia)', '(Papio)',
        '(Chlorocebus)', '(Nasua)', "Geoffroy's cat", '(elephant)', 'Ocelot', '(mollusk)',
        'camel', '(Vulpes)', '(Peromyscus)', 'feral swine', 'Hedgehog', '(Orcinus)', 
        '(Chrysomya)', '(Giraffa)', '(Macaca)', '(Panthera)', '(Gorilla)', '(Pongo)',
        '(Dipodomys)', '(Bison)', 'marmoset', '(Puma)', '(macaque)', '(Ailuropoda)',
        '(Nycticeius)', '(Odocoileus)', '(Myotis)', '(Procyon)', '(rufus)', '(Cervus)',
        '(Coati)', '(Marmota)', '(latrans)', '(clam)', '(Sylvicapra)', '(Otariinae)', 
        '(Camelus)', '(Okapia)', 'Cetacea', 'Serpentes', 'Mephitidae', 'Dasypodidae',
        'Sirenia', 'Phocidae', 'primate', 'primates', '(sika)', '(Ammotragus)', 
        '(Oreotragus)', '(Gazella)', '(gazella)', 'Cephalophus', '(Cephalophus)', 'moulfon', '(urial)', '(ibex)',
        '(goral)', '(slender-horned)', '(ibex)', '(roe)', 'bontebok', 'gerenuk', 'okapi',
        '(Nanger)', 'duiker', 'otter', 'orangutan', 'rhinoceros', 'Rhinoceros', '(Rucervus)', 
        '(Pseudois)', 'muntjac', '(wapiti)', '(Crotalus)', '(Antilocapra)', '(Hippotragus)',
        '(Rucervus)', '(Tragelaphus)', 'leopard', '(Eidolon)', '(Tursiops)', '(Desmodus rotundus)', 
        'Desmodus rotundus', 'raccoon', 'skunk', 'opossum', 'fox', 'wild boar', '(Delphinapterus)',
        'Paguma larvata', 'Canis lupus', '(bat)', '(gazelle)', 'Tadarida brasiliensis', 'Ovis canadensis',
        'Tayassu pecari', 'Tamandua mexicana', 'Buffalo', '(Buffalo)', 
        ]

ani_regex =  ('|'.join(ani))

# agricultural or domestic
agri = ['(cattle)' '(Bos)', '(taurus)', '(cow)', '(cow)', 'Cow', '(Cow)', '(bovine)', '(Bovine)',
        '(Calf)', '(cows)', '(Gallus)', 'gallus', '(Chicken)', 
        '(hen)', '(poultry)', 'Poultry', 'scrofa', '(porcine)', '(Pig)', '(pig)', '(swine)',
        '(Swine)', '(sheep)', 'aries', 'Canis lupus familiaris', '(Canis lupus familiaris)' ,'(ovine)', '(goat)', '(Dog)',
        '(dog)' '(Canine)', 'Canine', 'canine', '(familiaris)', 'familiaris' '(Equus)', '(equine)', '(horse)',
        '(duck)', '(Duck)', '(platyrhynchos)', '(goose)', '(cat)', '(catus)', '(turkey)',
        '(gallopavo)', '(rabbit)', '(rat)', '(Rat)', 'Cattle', 'Goat', 'Turkey',
        '(chicken)', 'chicken', '(chick)', 'chick', '(calf)', 'calf', 
        'dog', 'Bos', '(domestica)', 'domestica', 'domesticus', '(domesticus)',
        '(Ovine)', '(Pork)', '(Companion)', '(farm)', '(musculus)', 'donkey', 'steer', '(Cavia)', 
        '(Melopsittacus)', '(Vicugna)', '(Mustela)', '(Oryctolagus)', '(Capra)',
        '(Phasianus)', '(Chinchilla)', 'Meleagris', 'cockatoo', 'parrot', 'quail',
        'tortoise', 'llama', 'mink', 'Mink', 'psittacine', 'pigeon', 'mouse', '(Neovison)',
        '(Psittacidae)', 'Psittacidae', 'Murine', 'feline', 'Cricetinae', '(Neogale vison)',
        'Neogale vison']

agri_regex =  ('|'.join(agri))

# wildbirds
bird = ['(Oystercatcher)', '(Tern)', '(Sanderling)', '(Calidris)', '(Sandpiper)', 
        '(Plover)', '(Godwit)', '(Chroicocephalus)', '(Migratory)', '(Vultur)',
        '(Threskiornis)', '(Larus)', '( Spheniscus)', '(hornbill)', '(Swan)', '(swan)',
        'crow', 'Crow', '(Hirundo)', 'Seagull', 'wild bird', 'macaw', '(Rynchops)',
        'egret', '(Strigops)', '(Branta)', '(Haliaeetus)', '(Cygnus)', 'Cygnus',
        '(Quiscalus)', '(Catherpes)', '(Lagopus)', '(Passer)', '(Rupornis)',
        '(vulture)', '(owl)', '(ararauna)', 'Laridae', '(Polysticta)', '(Somateria)',
        '(Zenaida)', 'Pelecanus', '(Struthio)', '(Kea sp.)', 'Accipitridae', 'Anserini',
       '(Microcarbo)', '(Streptopelia)', 'swallow', 'frigatebird', 'Bubo', '(Ptyonoprogne)',
       '(Pica)', '(Corvus)', 'swallow']
        
bird_regex =  ('|'.join(bird))

missing = ['missing', '(missing)', 'Missing', '(Missing)', '(NOT)', '(Not)', '(not)', 'not',
           'Unknown', '(Unknown)', 'unknown', '(unknown)', 'none']

missing_regex =  ('|'.join(missing))

# sort dataframe using regex

sort = bs_raw.copy()

sort.loc[sort['host'].astype(str).str.contains(pat='r"' + human_regex + '"'), 'category'] = 'human'

sort.loc[sort['host'].astype(str).str.contains(pat='r"' + ani_regex + '"'), 'category'] = 'wild_animal'

sort.loc[sort['host'].astype(str).str.contains(pat='r"' + agri_regex + '"'), 'category'] = 'agricultural/domestic'

sort.loc[sort['host'].astype(str).str.contains(pat='r"' + bird_regex + '"'), 'category'] = 'wild_bird'

sort.loc[sort['host'].astype(str).str.contains(pat='r"' + missing_regex + '"'), 'category'] = 'missing'

sort['category'] = sort['category'].replace(np.nan, "other", regex=True)

sort.to_csv('../../results/host_sort_EC.csv')
