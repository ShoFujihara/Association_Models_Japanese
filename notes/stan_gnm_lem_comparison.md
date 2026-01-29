# Stan・gnm・LEM によるUnidiffモデルの比較

## 概要

Xie (1992) のFI_xモデル（Cross-nationally log-multiplicative full two-way interaction）を3つの手法で推定し、結果を比較した。

## データ

Yamaguchi (1987) / Xie (1992) の3カ国職業移動表データ
- US: N = 18,186
- Britain: N = 5,889
- Japan: N = 2,654
- 5×5×3 = 75セル

## モデル

FI_xモデル：
```
log μ_ijk = λ + λ_i^O + λ_j^D + λ_k^C + λ_ik^OC + λ_jk^DC + λ_ij^Diag×C + ψ_k × φ_ij
```

- O: Origin（父職業）, D: Destination（本人職業）, C: Country（国）
- ψ_k: 国ごとのスケーリングパラメータ（US=1に正規化）
- φ_ij: OD間の連関パターン（国間で共通）
- Diag×C: 対角セル（滞留者）を国ごとに飽和

## 結果比較

### スケーリングパラメータ（ψ）

| 手法 | US | Britain | Japan |
|------|-----|---------|-------|
| gnm | 1.0000 | 1.0398 | 0.7991 |
| LEM | 1.0000 | 1.0398 | 0.7989 |
| Stan MLE | 1.0000 | 1.1197 | 0.9052 |
| Stan MCMC | ~1.00 | ~1.10 | ~0.89 |

### モデル適合度（L²）

| 手法 | L² | df | p値 |
|------|-----|-----|-----|
| gnm | 30.94 | 20 | .056 |
| LEM | 30.94 | 20 | .056 |
| Stan MLE | 36.19 | - | - |

### 期待度数の相関

- gnm vs LEM: r = 1.0000（完全一致）
- gnm vs Stan MLE: r = 0.9999

## 実装詳細

### gnm

```r
gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C), O*D),
    family = poisson, data = d_xie)
```

### LEM

```
man 3
dim 3 5 5
mod {AB AC spe(BC,5a,A,c,2) spe(BC,1a,A,b)}
des [0 1 2]
```

- `spe(BC,5a,A,c,2)`: 対角セルのfull interaction（国ごと）
- `spe(BC,1a,A,b)`: Unidiff（乗法的層効果）

### Stan

cmdstanrの`$optimize()`でL-BFGS法による最尤推定。
モデル仕様は14-Bayes.qmd参照。

## 考察

### gnmとLEMの一致

gnmとLEMは完全に一致した結果を与える。これは両手法が同一のモデルを推定していることを示す。

### Stanで一致させるための注意点

**重要**: Stanでgnm/LEMと同じ結果を得るには、以下の点に注意が必要。

#### 1. 対角セルの完全飽和

gnmの`Diag*C`は15個の対角セル（5対角×3国）を完全に飽和させる。Stanでは：

```stan
// 誤: 階乗構造（パラメータ不足）
vector[D_max-1] alpha_diag_raw;           // 4パラメータ
matrix[D_max-1, K-1] alpha_diag_layer_raw; // 8パラメータ
// → 合計12パラメータ、15セルに不足

// 正: 完全飽和
vector[N_diag-1] alpha_diag_raw;  // 14パラメータ（N_diag=15）
```

#### 2. Unidiffは非対角セルのみに適用

対角セルにはUnidiff項を適用しない：

```stan
if (diag_cell_idx[n] > 0) {
  // 対角セル: 飽和パラメータのみ
  log_mu[n] += alpha_diag[diag_cell_idx[n]];
} else {
  // 非対角セル: Unidiff適用
  log_mu[n] += beta[layer_idx[n]] * psi[cell_idx[n]];
}
```

#### 3. 検証方法

モデル実装後、必ず以下を確認：
- 対角セルの期待度数が観測値と一致（飽和の確認）
- L²がgnm/LEMと一致
- スケーリングパラメータがgnm/LEMと一致

### 修正後の結果

| 手法 | US | Britain | Japan | L² |
|------|-----|---------|-------|-----|
| gnm | 1.0000 | 1.0398 | 0.7991 | 30.94 |
| Stan MLE（修正後） | 1.0000 | 1.0399 | 0.7991 | 31.00 |
| LEM | 1.0000 | 1.0398 | 0.7989 | 30.94 |

期待度数の相関: r = 1.0000（完全一致）

### ベイズ推定（MCMC）

修正後のモデルでベイズ推定も正しく動作する：

| パラメータ | gnm (MLE) | Stan (事後平均) | Stan (SD) |
|-----------|-----------|----------------|-----------|
| US | 1.000 | 1.000 | 0.000 |
| Britain | 1.040 | 1.030 | 0.071 |
| Japan | 0.799 | 0.795 | 0.087 |

- 事後平均はMLEとほぼ一致
- 信用区間が得られる（Britain: 0.91-1.15, Japan: 0.65-0.94）
- 期待度数の相関: r = 0.9999993

## ファイル

- MLE比較スクリプト: `scripts/compare_stan_mle_fixed.R`
- ベイズ推定テスト: `scripts/test_unidiff_bayes.R`
- LEM入力ファイル: `/Users/sf/Downloads/-old/lemfiles/examples/xie92_unidiff.INP`
- gnm実装: `11-Replication.qmd` (line 529)
- Stan実装（修正版）: `14-Bayes.qmd`

## 参考文献

- Xie, Y. (1992). The log-multiplicative layer effect model for comparing mobility tables. *American Sociological Review*, 57(3), 380-395.
- Yamaguchi, K. (1987). Models for comparing mobility tables: Toward parsimony and substance. *American Sociological Review*, 52(4), 482-494.
- Erikson, R., & Goldthorpe, J. H. (1992). *The Constant Flux*. Oxford: Clarendon Press.
