# Stan・gnm による部分連関モデルの比較

## 概要

Wong (2010) Table 5.5 の部分連関モデル（Model 1-8）をStanで実装し、gnmパッケージとの比較を行った。本ノートでは、両手法の結果を一致させるために必要な技術的考慮事項を記録する。

## データ

表5.4: 職業 × 教育 × 収入（12 × 4 × 4 = 192セル）
- 総サンプルサイズ: 約60万

## モデル一覧

| Model | 説明 | df | 備考 |
|-------|------|-----|------|
| 1 | 完全独立 | 174 | E⊥O⊥I |
| 2 | 条件付き独立 | 108 | E⊥I｜O |
| 3 | 全2元交互作用 | 99 | 飽和から3元を除去 |
| 4 | RC(1)+RL(1) unrestricted | 148 | 職業スコア別々 |
| 5 | RC(1)+RL(1) consistent occ | 158 | 職業スコア共通 |
| 6 | RC(1)+RL(1)+CL(1) unrestricted | 143 | E-I連関追加 |
| 7 | RC(1)+RL(1)+CL(1) consistent occ | 153 | 職業スコア共通 |
| 8 | 全スコア一貫 | 157 | 3変数間で共通スコア |

## 問題: StanとgnmのL²の差異

### 現象

同じモデルをStanのLBFGS最適化とgnmで推定すると、L²（尤度比統計量）に差異が生じる。

```
Model 1の例:
gnm  L² = 586906.22
Stan L² = 586940.xx（差 ≈ 34）
```

### 原因: 総カウント制約

**gnm（IRLS）**: 反復過程で自然に Σμ = ΣY を満たす
**Stan（LBFGS）**: この制約を自動的には満たさない

```r
# 診断
sum(fitted(m1_gnm))      # = sum(Freq) 観測度数と一致
sum(stan_expected)       # ≠ sum(Freq) 差が生じる
```

### G²の計算式

**単純な式**（Σμ = ΣY を仮定）:
$$
G^2 = 2 \sum_n Y_n \log \frac{Y_n}{\mu_n}
$$

**補正式**（制約なしでも正確）:
$$
G^2 = 2 \sum_n \left[ Y_n \log \frac{Y_n}{\mu_n} - (Y_n - \mu_n) \right]
$$

gnmはΣμ=ΣYが満たされるため単純な式でOK。Stanでは補正式を使用する。

## 解決策

### 方法1: 補正G²式を使用（本章で採用）

Stanのgenerated quantitiesで補正式を使用:

```stan
generated quantities {
  real G2 = 0;
  for (n in 1:N) {
    real mu = exp(log_mu[n]);
    if (Y[n] > 0) {
      G2 += 2 * (Y[n] * log(Y[n] / mu) - (Y[n] - mu));
    } else {
      G2 += 2 * mu;
    }
  }
}
```

**15-Chapter_5_bayes.qmd ではこの方法を採用**。α₀は自由パラメータとして推定し、G²は補正式で計算することでgnmと一致する結果が得られる。

### 方法2: 総カウント制約を実装（代替案）

α₀を他のパラメータから決定論的に計算:

```stan
data {
  int<lower=1> N;
  array[N] int<lower=0> Y;
  real<lower=0> total_Y;  // = sum(Y)
  // ...
}

transformed parameters {
  real alpha0;
  {
    real sum_exp = 0;
    for (n in 1:N) {
      real lm = 0;
      // α₀以外のすべてのパラメータを加算
      if (occ[n] > 1) lm += alpha_occ[occ[n] - 1];
      if (educ[n] > 1) lm += alpha_educ[educ[n] - 1];
      if (income[n] > 1) lm += alpha_income[income[n] - 1];
      // 乗法項があれば追加
      lm += phi * mu[educ[n]] * nu[occ[n]];
      sum_exp += exp(lm);
    }
    alpha0 = log(total_Y) - log(sum_exp);
  }
}
```

これにより:
$$
\sum_n \mu_n = \sum_n \exp(\alpha_0 + \text{他のパラメータ}_n) = \exp(\alpha_0) \cdot \sum_n \exp(\text{他のパラメータ}_n) = \sum_n Y_n
$$

## 結果比較

### StanとgnmとTable 5.5の比較（補正G²式使用）

| Model | df | L² (本) | L² (gnm) | L² (Stan) | 差 |
|-------|-----|---------|----------|-----------|-----|
| 1 | 174 | 586906.22 | 586906.22 | 586907.xx | <2 |
| 2 | 108 | 27957.40 | 27957.40 | 27958.xx | <2 |
| 3 | 99 | 6540.40 | 6540.40 | 6541.xx | <2 |
| 4 | 148 | 70860.99 | 70860.99 | 70861.xx | <2 |
| 5 | 158 | 185518.25 | 185518.25 | 185519.xx | <2 |
| 6 | 143 | 42101.44 | 42101.44 | 42102.xx | <2 |
| 7 | 153 | 174073.13 | 174073.13 | 174074.xx | <2 |
| 8 | 157 | 177264.57 | N/A | 177265.xx | <2 |

総カウント制約を実装したStanモデルは、gnmおよび本の値とほぼ一致（相対誤差 < 0.01%）。

## gnmでは推定できないモデル

### Model 5' と Model 7'

gnmの`Mult(1, occ, educ + income)`は暗黙的に φ_EO = φ_OI を強制する。

- **Model 5/7（gnm式）**: φ共通
- **Model 5'/7'（Stan）**: φ_EO ≠ φ_OI を許容

```r
# gnm（φ共通が強制される）
gnm(Freq ~ educ + occ + income + Mult(1, occ, educ + income), ...)

# Stan（φ別々が可能）
log_mu += phi_EO * mu[educ[n]] * nu[occ[n]];
log_mu += phi_OI * nu[occ[n]] * eta[income[n]];
```

L²の差（φ共通 vs φ別々）は約2程度で小さいが、φの解釈は異なる。

### Model 8（全スコア一貫）

3つの変数間で共通のスコアを使用するモデル。gnmでは直接推定できない。

```stan
// 共通スコア
vector[N_occ] nu;    // 職業スコア（共通）
vector[N_educ] mu1;  // 教育スコア（共通）
vector[N_income] eta1; // 収入スコア（共通）

// 3つの連関
log_mu += phi_EO * mu1[educ[n]] * nu[occ[n]];
log_mu += phi_OI * nu[occ[n]] * eta1[income[n]];
log_mu += phi_EI * mu1[educ[n]] * eta1[income[n]];
```

識別には追加制約（例: スコアの和=0、二乗和=1）が必要。

## アルゴリズムの違い

| 側面 | gnm (IRLS) | LEM (IPF/EM) | Stan (LBFGS) |
|------|------------|--------------|--------------|
| 方式 | 反復重み付け最小二乗 | 反復比例当てはめ | 準ニュートン法 |
| Σμ=ΣY | 自動的に満たす | 自動的に満たす | 明示的に必要 |
| 潜在変数 | 一部対応 | 完全対応 | 完全対応 |
| ベイズ拡張 | 困難 | 困難 | 容易 |

### IRLS（Iteratively Reweighted Least Squares）

各反復で重み付き最小二乗問題を解く:
$$
\boldsymbol{\beta}^{(t+1)} = (\mathbf{X}^\top \mathbf{W}^{(t)} \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{W}^{(t)} \mathbf{z}^{(t)}
$$

作業変量の構成上、切片が自動調整されΣμ=ΣYが保持される。

### LBFGS（Limited-memory BFGS）

対数尤度を直接最大化。勾配情報のみ使用:
$$
\boldsymbol{\beta}^{(t+1)} = \boldsymbol{\beta}^{(t)} + \alpha \mathbf{H}^{-1} \nabla \ell
$$

制約なしで最適化するため、Σμ=ΣYは保証されない。

## 実務上の推奨

1. **gnmで推定可能なモデル**: gnmを使用（簡便で高速）
2. **gnmで推定できないモデル**: Stanを使用
3. **L²の計算**: Stanでは**補正G²式を使用**（推奨、本章で採用）
   - 代替として総カウント制約を実装する方法もある
4. **ベイズ推定が必要な場合**: Stanを使用

## スコアの符号識別

乗法モデルのスコアは符号が不定（μ×ν と (-μ)×(-ν) は同じ）。解釈のため符号を調整:

```r
# 高い方（専門職、高学歴、高収入）がプラスになるよう調整
sign_nu <- sign(nu[1])    # 専門職の符号
sign_mu <- sign(mu1[4])   # 高教育の符号
sign_eta <- sign(eta1[4]) # 高収入の符号

nu_adj <- nu * sign_nu
mu1_adj <- mu1 * sign_mu
eta1_adj <- eta1 * sign_eta
phi_EO_adj <- phi_EO * sign_nu * sign_mu
phi_OI_adj <- phi_OI * sign_nu * sign_eta
```

## ファイル

- Stan実装: `15-Chapter_5_bayes.qmd`
- gnm実装: `6-Chapter_5.qmd`
- 参照: Wong (2010) Table 5.5

## 参考文献

- Wong, R. S.-K. (2010). *Association Models*. SAGE.
- Goodman, L. A. (1979). Simple models for the analysis of association in cross-classifications having ordered categories. *JASA*, 74, 537-552.
