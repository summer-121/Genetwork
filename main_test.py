from lib_back import Pub_Analysis
from lib_back import Slope
from lib_back import HelpGPT

gene = "CREB"

if __name__ == "__main__":
    analyzer = Pub_Analysis(email="1018jjkk@gmail.com")
    yearly_data = analyzer.get_yearly_publications(gene, 2011, 2020)

    slope = Slope(yearly_data)
    df_with_change = slope.compute_rate_of_change()
    significant_years = slope.detect_significant_changes(threshold=2.0)

    gpt_helper = HelpGPT(model="gpt-5")
    result_df = gpt_helper.analyze_changes(df_with_change, significant_years, gene)

    print(result_df)
