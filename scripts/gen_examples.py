"""Generate example data files for all 4 model types."""
import pandas as pd
from pathlib import Path

BASE = Path(__file__).parent.parent
EXAMPLES = BASE / "data" / "examples"
EXAMPLES.mkdir(parents=True, exist_ok=True)

# MDI: top diseases
mdi = pd.read_csv(BASE / "data/databases/mdi_database.csv",
                   usecols=["Disease_Name", "Disease_ID"])
top_d = (mdi.groupby(["Disease_Name", "Disease_ID"]).size()
         .reset_index(name="count").nlargest(10, "count"))
print("=== MDI top diseases ===")
print(top_d.to_string())

# MMI: top microbes
mmi = pd.read_csv(BASE / "data/databases/mmi_database.csv",
                   usecols=["Microbe_Name", "Taxonomy_ID"])
top_m = (mmi.groupby(["Microbe_Name", "Taxonomy_ID"]).size()
         .reset_index(name="count").nlargest(10, "count"))
print("\n=== MMI top microbes ===")
print(top_m.to_string())

# MDrI: top drugs
mdri = pd.read_csv(BASE / "data/databases/mdri_database.csv",
                    usecols=["Drug_Name", "DrugBank_ID"])
top_dr = (mdri.groupby(["Drug_Name", "DrugBank_ID"]).size()
          .reset_index(name="count").nlargest(10, "count"))
print("\n=== MDrI top drugs ===")
print(top_dr.to_string())

# Also get example metabolites for each model type — common metabolites
met_mdi = pd.read_csv(BASE / "data/databases/mdi_database.csv",
                       usecols=["Metabolite_Name", "HMDB_ID"])
met_mdi_top = (met_mdi.groupby(["Metabolite_Name", "HMDB_ID"]).size()
               .reset_index(name="count").nlargest(10, "count"))
print("\n=== MDI top metabolites ===")
print(met_mdi_top.to_string())

met_mmi = pd.read_csv(BASE / "data/databases/mmi_database.csv",
                       usecols=["Metabolite_Name", "HMDB_ID"])
met_mmi_top = (met_mmi.dropna(subset=["Metabolite_Name"])
               .groupby(["Metabolite_Name"]).size()
               .reset_index(name="count").nlargest(10, "count"))
print("\n=== MMI top metabolites ===")
print(met_mmi_top.to_string())

met_mdri = pd.read_csv(BASE / "data/databases/mdri_database.csv",
                        usecols=["Metabolite_Name", "HMDB_ID"])
met_mdri_top = (met_mdri.groupby(["Metabolite_Name", "HMDB_ID"]).size()
                .reset_index(name="count").nlargest(10, "count"))
print("\n=== MDrI top metabolites ===")
print(met_mdri_top.to_string())
