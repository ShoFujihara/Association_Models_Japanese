# 12-Review.qmd 参考文献の詳細

本ファイルは12-Review.qmd（研究の動向）で引用している論文の正確な情報をまとめたものである。

## 発見された問題点

### 重大な誤り

1. **@breen2010** - シミュレーション研究ではなく、「Educational Expansion and Social Mobility in the 20th Century」（教育拡大と社会移動）
2. **@pfeffer2015** - 疎データの研究ではなく、「How Has Educational Expansion Shaped Social Mobility Trends in the United States」（教育拡大と社会移動のトレンド）
3. **@sobel1985** - bibliographyに2つの論文が同じキーで存在。SHDモデルの論文は「Exchange, Structure, and Symmetry in Occupational Mobility」(AJS)であり、「Social Mobility and Fertility Revisited」(ASR)ではない
4. **@hout1984** - bibliographyの論文が「Occupational Mobility of Black Men」になっているが、正しくは「Status, Autonomy, and Training in Occupational Mobility」

### 説明の修正が必要な箇所

5. **@zhou2015** - 「連関係数の性質を比較」ではなく、「Shrinkage Estimation of Log-odds Ratios for Comparing Mobility Tables」（対数オッズ比の縮小推定）
6. **@karlson2023a** - 「オッズ比の解釈」一般ではなく、「Marginal Odds Ratios」（周辺オッズ比）に特化
7. **@mize2019a** - 「限界効果」ではなく、「Nonlinear Interaction Effects」（非線形交互作用効果）

---

## 論文詳細

### 連関モデルの基礎

#### @goodman1979a
- **タイトル**: Simple Models for the Analysis of Association in Cross-Classifications Having Ordered Categories
- **著者**: Leo A. Goodman
- **雑誌**: Journal of the American Statistical Association, Vol. 74, pp. 537-552, 1979
- **内容**: 順序カテゴリを持つクロス表における連関分析のための単純モデルを提案。隣接する行と列からなる2×2小表のオッズ比として連関を測定。RC連関モデルの基礎論文。
- **URL**: https://www.tandfonline.com/doi/abs/10.1080/01621459.1979.10481650

#### @duncan1979
- **タイトル**: How Destination Depends on Origin in the Occupational Mobility Table
- **著者**: Otis Dudley Duncan
- **雑誌**: American Journal of Sociology, Vol. 84, pp. 793-804, 1979
- **内容**: 1949年イギリスの職業移動表を再分析し、「達成」の観点を強調する2つのモデルを提示。数学的簡潔さ、統計的適合度、概念的明確さ、比較研究への適用性において魅力的なモデルを提案。
- **URL**: https://www.journals.uchicago.edu/doi/10.1086/226861

### SATモデル

#### @hout1984（要修正）
- **正しいタイトル**: Status, Autonomy, and Training in Occupational Mobility
- **著者**: Michael Hout
- **雑誌**: American Journal of Sociology, Vol. 89, pp. 1379-1409, 1984
- **内容**: 社会経済的背景が職業達成に与える効果は確立されているが、地位の効果だけでは出身と到達の連関を説明しきれない。労働者に与えられる**自律性（autonomy）**と必要とされる**訓練の専門性（training）**も移動に重要。自律性と訓練は特に不動（immobility）に重要。OCG 1962年・1973年データに適用。
- **URL**: https://www.journals.uchicago.edu/doi/abs/10.1086/228020
- **注意**: bibliographyの@hout1984が「Occupational Mobility of Black Men」を指している場合は修正が必要

#### @hout1988
- **タイトル**: More Universalism, Less Structural Mobility: The American Occupational Structure in the 1980s
- **著者**: Michael Hout
- **雑誌**: American Journal of Sociology, Vol. 93, pp. 1358-1400, 1988
- **内容**: 1972年から1985年にかけて、男女の社会経済的出身と到達の連関が3分の1減少。このトレンドは大卒労働者の増加と関連。大卒は背景地位の効果を相殺するため、大卒が増えるほど出身-到達連関は弱まる。構造移動の減少が階層構造の開放性増加を相殺し、全体の移動率は不変。
- **URL**: https://www.journals.uchicago.edu/doi/10.1086/228904

### SHDモデル

#### @sobel1985（正しい論文）
- **タイトル**: Exchange, Structure, and Symmetry in Occupational Mobility
- **著者**: Michael E. Sobel, Michael Hout, Otis Dudley Duncan
- **雑誌**: American Journal of Sociology, Vol. 91, pp. 359-372, 1985
- **内容**: 従来の交換移動と構造移動の概念化が対数線形モデルのパラメータと適切に対応していなかった問題を指摘。準対称モデル（QS）のパラメータと概念の対応関係を確立。
  - **交換移動（reciprocated mobility）**: 職業カテゴリ間の等しいフロー（対称的な部分）
  - **構造移動**: 周辺分布の異質性が出身に一様に作用する効果
  - QSまたはその特殊ケースが成立する場合に、パラメータと概念の対応関係が成立
  - アドホックでない、パラメトリックな構造移動指標の構築が可能に
- **URL**: https://www.journals.uchicago.edu/doi/abs/10.1086/228281
- **注意**: bibliographyに「Social Mobility and Fertility Revisited」(ASR)も同じキーで存在する可能性あり

### 対数乗法層効果モデル・Unidiff

#### @yamaguchi1987
- **タイトル**: Models for Comparing Mobility Tables: Toward Parsimony and Substance
- **著者**: Kazuo Yamaguchi
- **雑誌**: American Sociological Review, Vol. 52, pp. 482-494, 1987
- **内容**: 一様層効果モデル（uniform layer effect model）を提案。2つの移動表間の出身-到達連関の差を1パラメータで記述。「垂直移動」の差異を分析するための1パラメータ検定を提供。
- **URL**: https://www.jstor.org/stable/pdf/2095293.pdf

#### @wong1990
- **タイトル**: Understanding Cross-National Variation in Occupational Mobility
- **著者**: Raymond Sin-Kwok Wong
- **雑誌**: American Sociological Review, Vol. 55, pp. 560-573, 1990
- **内容**: 国際比較移動研究の基礎的論文。対数乗法層効果モデル（Unidiff）を用いて国間の職業移動パターンの違いを分析。

#### @erikson1992
- **タイトル**: The Constant Flux: A Study of Class Mobility in Industrial Societies
- **著者**: Robert Erikson, John H. Goldthorpe
- **出版社**: Oxford: Clarendon Press, 1992
- **内容**: 近代産業社会における社会移動の研究。第二次世界大戦後の「長期ブーム」期における欧州諸国（西欧・東欧）の経験を中心に、米国・豪州・日本の章も含む。CmSFモデル、コアフルイディティ、FJH仮説などを扱う。Unidiffモデルの適用事例。
- **URL**: https://books.google.com/books/about/The_constant_flux.html?id=WfYDAQAAIAAJ

#### @xie1992
- **タイトル**: The Log-Multiplicative Layer Effect Model for Comparing Mobility Tables
- **著者**: Yu Xie
- **雑誌**: American Sociological Review, Vol. 57, pp. 380-395, 1992
- **内容**: 対数乗法層効果モデルを提案。表間の出身-到達連関の変動を、共通の連関パターンと表固有のパラメータの対数乗法積として制約。Yamaguchi (1987)の一様層効果モデルより柔軟。比較研究での簡潔さと解釈可能性が魅力。
- **URL**: https://www.jstor.org/stable/2096242

### 回帰型モデル

#### @goodman1998
- **タイトル**: Statistical Methods and Graphical Displays for Analyzing How the Association between Two Qualitative Variables Differs among Countries, among Groups, or over Time: A Modified Regression-Type Approach
- **著者**: Leo A. Goodman, Michael Hout (with comments by Yu Xie, Kazuo Yamaguchi)
- **雑誌**: Sociological Methodology, Vol. 28, pp. 175-230, 1998
- **内容**: 2つの質的変数間の連関が、国・グループ・時間によってどのように異なるかを分析するための統計的方法とグラフ表示。修正された回帰型アプローチを提案。
- **URL**: https://www.jstor.org/stable/270967

#### @goodman2001
- **タイトル**: Statistical Methods and Graphical Displays for Analyzing How the Association between Two Qualitative Variables Differs among Countries, among Groups, or over Time: Part II: Some Exploratory Techniques, Simple Models, and Simple Examples
- **著者**: Leo A. Goodman, Michael Hout
- **雑誌**: Sociological Methodology, Vol. 31, pp. 189-221, 2001
- **内容**: 1998年論文の続編。探索的技法、単純モデル、単純な例を提示。

### 連関係数・要約統計量

#### @zhou2015
- **タイトル**: Shrinkage Estimation of Log-odds Ratios for Comparing Mobility Tables
- **著者**: Xiang Zhou
- **雑誌**: Sociological Methodology, Vol. 45, 2015
- **内容**: 移動表比較のための対数オッズ比の縮小推定量を提案。経験ベイズの枠組みで、複数の表から「強さを借りる」ことで推定効率を改善。数値シミュレーションでMLEより優れた性能を示す。Altham指数の調整推定量も構築し、社会流動性の全体的程度を測る際にUnidiffモデルの結果とより一致。
- **URL**: https://journals.sagepub.com/doi/10.1177/0081175015570097
- **注意**: 「連関係数の性質を比較」という記述は不正確。縮小推定に特化。

#### @bouchet-valat2022
- **タイトル**: General Marginal-free Association Indices for Contingency Tables: From the Altham Index to the Intrinsic Association Coefficient
- **著者**: Milan Bouchet-Valat
- **雑誌**: Sociological Methods & Research, Vol. 51(1), pp. 203-236, 2022
- **内容**: 周辺自由な連関指数の一般的枠組み。オッズ比、Altham指数、内在的連関係数、対数乗法モデルの係数の直接的関係を明らかにする。0から1の間で変動する正規化版を考案し、相関係数に似た解釈を提供。149の欧州地域における教育・社会経済的同類婚の事例で例示。
- **URL**: https://journals.sagepub.com/doi/10.1177/0049124119852389

### 因果推論との関係

#### @yamaguchi2012
- **タイトル**: Log-linear Causal Analysis of Cross-classified Categorical Data
- **著者**: Kazuo Yamaguchi
- **雑誌**: Sociological Methodology, Vol. 42, 2012
- **内容**: 対数線形・対数乗法の連関分析の主な限界は因果分析との関連の欠如であった。セミパラメトリック回帰モデルを仮定し、調整済みクロス表においてXのYへの因果効果を保持する新しい方法を導入。交絡変数Vの効果を特定しないセミパラメトリックロジット・多項ロジット回帰モデルの適用方法も紹介。
- **URL**: https://journals.sagepub.com/doi/10.1177/0081175012460661

#### @kuha2010
- **タイトル**: Path Analysis for Discrete Variables: The Role of Education in Social Mobility
- **著者**: Jouni Kuha, John H. Goldthorpe
- **雑誌**: Journal of the Royal Statistical Society. Series A, Vol. 173(2), pp. 351-369, 2010
- **内容**: 教育達成が世代間社会移動に与える効果の大きさを評価。一部の変数がカテゴリカルな場合でも直接効果と間接効果を推定できる一般的パス分析法を提案。効果を平均差で表現した場合に正確な加法分解を提供し、対数オッズ比では近似的だが通常かなり正確。英国調査データの社会移動分析で例示。
- **URL**: https://eprints.lse.ac.uk/31253/

### オッズ比・限界効果

#### @karlson2023a
- **タイトル**: Marginal Odds Ratios: What They Are, How to Compute Them, and Why Sociologists Might Want to Use Them
- **著者**: Kristian Bernt Karlson, Ben Jann
- **雑誌**: Sociological Science, Vol. 10, pp. 332-347, 2023
- **内容**: 社会学者がオッズ比から離れ、平均限界効果を報告することが増えている。周辺オッズ比を導入し、オッズ比の使用を復活させることを目指す。従来のオッズ比と異なり、周辺オッズ比は省略された共変量の影響を任意に受けない。潜在的アウトカムで定義し、平均限界効果との密接な関係、従来のオッズ比に対する利点を議論。
- **URL**: https://sociologicalscience.com/articles-v10-10-332/
- **注意**: 「オッズ比の解釈」一般ではなく、**周辺オッズ比**に特化した論文

#### @mize2019a
- **タイトル**: Best Practices for Estimating, Interpreting, and Presenting Nonlinear Interaction Effects
- **著者**: Trenton D. Mize
- **雑誌**: Sociological Science, Vol. 6, pp. 81-117, 2019
- **内容**: 非線形効果と交互作用効果の両方を含む分析について、推定・解釈・提示のベストプラクティスを提示。線形回帰における非線形効果と、カテゴリカルアウトカムモデル（二値ロジット/プロビット中心）における非線形効果を扱う。非線形交互作用効果を線形交互作用効果と同じように扱うことの深刻な結果について警告。
- **URL**: https://sociologicalscience.com/articles-v6-4-81/
- **注意**: 「限界効果」一般ではなく、**非線形交互作用効果**に特化した論文

### 対応分析との関係

#### @vanderheijden1985
- **タイトル**: Correspondence Analysis Used Complementary to Loglinear Analysis
- **著者**: Peter G. M. van der Heijden, Jan de Leeuw
- **雑誌**: Psychometrika, Vol. 50, pp. 429-447, 1985
- **内容**: 対数線形分析と対応分析がクロス表の分解に異なる方法を提供することを示す。対応分析が2つの行列（それぞれ特定の対数線形モデルに従う）の差の分解として見なせる場合があることを示す。これらの場合、対応分析の解は対数線形モデル間の差として解釈可能。Escofierによる対応分析の一般化も議論。
- **URL**: https://link.springer.com/article/10.1007/BF02296262

#### @vanderheijden1988
- **タイトル**: Comment on "Correspondence analysis used complementary to loglinear analysis"
- **著者**: Peter G. M. van der Heijden, Keith J. Worsley
- **雑誌**: Psychometrika, Vol. 53, pp. 287-291, 1988
- **内容**: 1985年論文へのコメント。

---

## 反実仮想的分解（Counterfactual Decomposition）

@breen2010 と @pfeffer2015 は**反実仮想的分解（counterfactual decomposition）**の手法を用いた研究である。これらは「シミュレーション研究」ではなく、教育拡大が社会移動に与える因果的効果を分解する方法論的貢献である。

### 分解手法の概要

Breen (2010) が開発し、Pfeffer & Hertel (2015) が発展させた分解手法は、観察された社会流動性のトレンドを、教育の役割に関する特定の仮定に基づく反実仮想的トレンドと比較する。

**3つのメカニズム**:
1. **構成効果（Compositional Effect, OED）**: 高学歴者が増加することで、出身階層の直接的影響が弱まる
2. **教育における階層不平等の変化（COE）**: コーホート間で教育達成の階層格差が変化する
3. **教育収益の変化（CED）**: コーホート間で教育の階層達成への効果が変化する

### @breen2010
- **タイトル**: Educational Expansion and Social Mobility in the 20th Century
- **著者**: Richard Breen
- **雑誌**: Social Forces, Vol. 89(2), pp. 365-388, 2010
- **内容**: 反実仮想的分解手法を開発。世代間社会階層移動における教育の役割を測定するツールを提示。英国・スウェーデン・ドイツの20世紀における教育拡大の効果を分析し、3カ国すべてで教育拡大がより大きな社会移動を促進したことを発見。
- **方法論的貢献**: 観察されたトレンドと反実仮想的トレンドを比較し、各メカニズムの相対的重要性を明らかにする分解手法
- **URL**: https://academic.oup.com/sf/article-abstract/89/2/365/2235230

### @pfeffer2015
- **タイトル**: How Has Educational Expansion Shaped Social Mobility Trends in the United States?
- **著者**: Fabian T. Pfeffer, Florian R. Hertel
- **雑誌**: Social Forces, Vol. 94(1), pp. 143-180, 2015
- **内容**: Breen (2010) の分解手法を用いて、1972-2012年のGSSデータ（N=14,588-14,608）を分析。社会階層移動の緩やかだが漸進的な増加は、ほぼ**構成効果**によるものと発見。教育分布は大きく変化したにもかかわらず、教育における階層不平等は安定し、教育収益にも一貫したトレンドがなかった。
- **URL**: https://pmc.ncbi.nlm.nih.gov/articles/PMC4543307/

---

## 修正アクション

### 完了済み（12-Review.qmd）

1. **@hout1984 → @hout1984a に変更**: bibliographyには両方のHout 1984論文が存在
   - @hout1984 = "Occupational Mobility of Black Men" (ASR)
   - @hout1984a = "Status, Autonomy, and Training in Occupational Mobility" (AJS) ← **SATモデルの論文はこちら**
2. **シミュレーション研究セクションを削除**: @breen2010と@pfeffer2015は教育拡大と社会移動の研究であり、シミュレーション研究ではない
3. **各論文の説明を正確な内容に修正**
4. **@breen2010と@pfeffer2015を「Unidiffモデルを用いた実証研究」として追加**: これらはUnidiffモデルを用いて教育拡大と社会移動の関係を分析した実証研究であるため、Unidiffセクションの下に新しいサブセクションとして追加

### 要確認（bibliography）

4. **@sobel1985のキー重複問題**: bibliographyに同じキー@sobel1985で2つの異なる論文が存在
   - "Social Mobility and Fertility Revisited" (ASR)
   - "Exchange, Structure, and Symmetry in Occupational Mobility" (AJS) ← **SHDモデルの論文はこちら**
   - 現状ではQuartoがどちらを参照するか不定。片方のキーを変更することを推奨（例: sobel1985 → sobel1985a）
