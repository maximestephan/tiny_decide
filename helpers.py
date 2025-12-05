from rdkit import Chem
import hashlib
from rdkit.Chem import Descriptors

def smiles_to_canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalid : {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)

POOL_SIZE = 26**2 * 10**6  # 676 000 000 possibilities

def canonical_string_to_id(canon: str, prefix: str = "MU") -> str:
    digest = hashlib.sha256(canon.encode("utf-8")).digest()
    n = int.from_bytes(digest, byteorder="big")
    n = n % POOL_SIZE

    num_part = n % 10**6          # 0 .. 999999
    letter_index = n // 10**6     # 0 .. 26^2 - 1

    first_letter_index = letter_index // 26
    second_letter_index = letter_index % 26

    first_letter = chr(ord("A") + first_letter_index)
    second_letter = chr(ord("A") + second_letter_index)

    return f"{prefix}{first_letter}{second_letter}{num_part:06d}"

def smiles_to_id(smiles: str, prefix: str = "MU"):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        print(f"SMILES invalid : {smiles}")
        return "MU0"
        #raise ValueError(f"SMILES invalide : {smiles}")
    canon = Chem.MolToSmiles(mol, canonical=True)
    return canonical_string_to_id(canon, prefix=prefix)

dictonnary = {}
dictonnary_smiles = {}
def test_id(input_file, output_file, error_file):
    with open(input_file, "r", encoding="utf-8") as fin, \
         open(output_file, "w", encoding="utf-8") as fout, \
         open(error_file, "w", encoding="utf-8") as ferr:
        for line in fin:
            line = line.rstrip("\n")
            column = line.split("\t")
            id, smiles, mol = smiles_to_id(column[1])
            if id == "MU0":
                pass
                #ferr.write(id + "\t" +  column[0] +   "\t"  +  column[1] + "\n")
            if id not in dictonnary:
                dictonnary[id] = column[0]
                dictonnary_smiles[id] = smiles
            else:
                if dictonnary_smiles[id] != smiles:
                    ferr.write(id + "\t" + dictonnary[id] + "\t" +  column[0]  + "\t" + dictonnary_smiles[id] + "\t" +  column[1] +"\n" )
            fout.write(id + "\t" +  column[0] + "\n")

#test_id("smiles.txt","id.txt","id_error.txt")
    
    
PESTICIDE_DESCRIPTOR_FUNCS = {
    "MolWt": Descriptors.MolWt,
    "CalcLogP": Descriptors.MolLogP,
    "TPSA": Descriptors.TPSA,  
    "NumHAcceptors": Descriptors.NumHAcceptors,
    "NumHDonors": Descriptors.NumHDonors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "FractionCSP3": Descriptors.FractionCSP3,
    "NumValenceElectrons": Descriptors.NumValenceElectrons,
    "NumHeavyAtoms": Descriptors.HeavyAtomCount,
    "RingCount": Descriptors.RingCount,
    "NumAromaticRings": Descriptors.NumAromaticRings,
    "NumAliphaticRings": Descriptors.NumAliphaticRings,
    "NumHeteroatoms": Descriptors.NumHeteroatoms,
}
HALOGEN_PATTERN = Chem.MolFromSmarts("[F,Cl,Br,I]")

def count_halogens(mol: Chem.Mol) -> int:
    if HALOGEN_PATTERN is None or mol is None:
        return 0
    return len(mol.GetSubstructMatches(HALOGEN_PATTERN))

def compute_pesticide_descriptors_from_mol(mol: Chem.Mol) -> dict:
    if mol is None:
        raise ValueError("mol est None")

    Chem.SanitizeMol(mol, catchErrors=True)

    descs = {}

    for name, func in PESTICIDE_DESCRIPTOR_FUNCS.items():
        try:
            descs[name] = func(mol)
        except Exception:
            descs[name] = None

    try:
        descs["FormalCharge"] = mol.GetFormalCharge()
    except Exception:
        descs["FormalCharge"] = None

    try:
        descs["NumHalogens"] = count_halogens(mol)
    except Exception:
        descs["NumHalogens"] = None

    return descs


def all_descriptors(mol):
    results = {}
    for name, func in Descriptors._descList:
        try:
            results[name] = func(mol)
        except Exception as e:
            results[name] = None
    return results

def chemistrify(input_file, output_file):
    with open(input_file, "r", encoding="utf-8") as fin, \
         open(output_file, "w", encoding="utf-8") as fout:
        fout.write("id" + "\t" +  "MUN" + "\t" + "smiles" )
        dummy_mol = Chem.MolFromSmiles("C")
        descs = compute_pesticide_descriptors_from_mol(dummy_mol)
        for k, v in descs.items():
            fout.write("\t" + k)
        fout.write("\n")
        for line in fin:
            line = line.rstrip("\n")
            column = line.split("\t")
            name =  column[0]
            smiles = column[1]
            mol = Chem.MolFromSmiles(smiles)
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            id = canonical_string_to_id(canonical_smiles)
            fout.write(id + "\t" + name + "\t" +  smiles )
            #logP =  round(Descriptors.MolLogP(mol),1)
            descs = compute_pesticide_descriptors_from_mol(mol)
            for k, v in descs.items():
                fout.write("\t" + str(v))
            fout.write("\n")
            
#chemistrify("./data/pesticides.smiles","./data/pesticides_chemistry.csv")

def round_to_one_significant_decimal(x):
    if "e" in str(f"{x:.3g}"):
        if x > 1 or x < -1:
            return int(x)
        else:
            decimals = -floor(log10(abs(x)))
            as_str = f'{{:.{decimals}f}}'.format(x)
            chunks = []
            start_chunk = as_str.find('.') + 4
            chunks.append(as_str[:start_chunk])
            for i in range(start_chunk, len(as_str), 3):
                chunks.append(as_str[i:i+3])
            return ' '.join(chunks)

    return f"{x:.3g}"
    #return round(x, -int(floor(log10(abs(x)))))
