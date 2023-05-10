# 『カテゴリカルデータの連関モデル』のサポートページ

Author: 藤原翔（東京大学社会科学研究所）
Email: sho.fujihara [atmark] iss.u-tokyo.ac.jp

ここは，Wong, Raymond Sin-Kwok. 2010. Association Models. Thousand Oaks: SAGE. の訳書である『カテゴリカルデータの連関モデル』（共立出版）のサポートページです．

# 目次
<!-- 目次部分(リンクになるところ) -->
1. [正誤表](#anchor1)
1. [Rの使い方](#anchor2)
1. [重要な対数線形・対数乗法モデルの再現](#anchor3)
1. [『カテゴリカルデータの連関モデル』](#anchor4)

<a id="anchor1"></a>
# 正誤表
正誤表については見つかり次第，ここで報告します．

また翻訳作業の中でWong (2010)で見つかった誤記についてもここでフォローしたいと思います．

- ζ*i. = 1 for all i -> ζ*i. = i for all

<a id="anchor2"></a>
# Rの使い方
- 基本的なRの使用方法を説明します．


## Rのインストール {-}

- CRAN（The Comprehensive R Archive Network）のページ（ https://cran.r-project.org ）
から最新版のRをダウンロードし，インストールてください．
- 必要であればRStudioもダウンロードし，インストールしてください．

```
a <- 1:5
```

# LEMの使い方


## $l_{\rm EM}$のインストール {-}

- $l_{\rm EM}$はJeroen K. Vermunt氏のホームページ（ https://jeroenvermunt.nl/ ）からダウンロードできる．
- Software, User Manuals, and VideoのVermunt, J.K  (1997). LEM 1.0: A general program for the analysis of categorical data. Tilburg: Tilburg University. (download) の「 (download) 」をクリック．
- Windowsで動く．

# 基本的なクロス表の分析方法

- 対数線形モデル，対数乗法モデル，また連関モデルを分析する前に，ここでかんたんなクロス表の分析方法について説明します．

<a id="anchor3"></a>
# 重要な対数線形・対数乗法モデルの再現
3つの重要な論文のモデルを再現します．
  - Duncan, Otis Dudley. 1979. “How Destination Depends on Origin in the Occupational Mobility Table.” *The American Journal of Sociology* 84(4):793–803. 
  - Yamaguchi, Kazuo. 1987. “Models for Comparing Mobility Tables: Toward Parsimony and Substance.” *American Sociological Review* 52(4):482-494.
  - Xie, Yu. 1992. “The Log-Multiplicative Layer Effect Model for Comparing Mobility Tables.” *American Sociological Review* 57(3):380-395.

Duncan (1979) については`uniform_association_model.R`で，Yamaguchi (1987)については，`uniform_association_model.R`で，そしてXie (1992)については，`uniform_layer_effect_model.R`で行います．


<a id="anchor4"></a>
# 『カテゴリカルデータの連関モデル』
サポートサイトとしては，Raymond Sin-Kwok WongによるStudent Study Site for Association Models （https://studysites.sagepub.com/wongstudy/ ）があります．
ここでは本書で十分に説明されなかった部分を中心にフォローしたいと思います．

以下のサイトにRのコードと結果をまとめています．
https://shofujihara.github.io/Association_Models/index.html

