import pandas as pd


if __name__ == "__main__":
    df = pd.read_csv("workflow/results/nullomer_statistics_summarized.csv", sep=',')
    df_json = pd.read_json("workflow/data/organisms.json").transpose()
    df_json.index = df_json.index.astype(str).str.replace('.','_', regex=False).str.strip()
    df_new = df.merge(df_json, left_on='organism', right_index=True, how='left')
    df_new['trivial_porc'] = (df_new['trivial']/df_new['counter'])*100
    print(df_new)
    df_new.to_csv("workflow/results/filtered_nullomer_statistics_summarized.csv", sep=',')
