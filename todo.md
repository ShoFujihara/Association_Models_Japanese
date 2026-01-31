# TODO - Association Models Japanese

## 現在のステータス

最終更新: 2026-01-31

## 完了済み

### 各章の校正・修正
- [x] 1-Introduction.qmd - typo修正、一様連関モデルのコード修正
- [x] 2-Chapter_1.qmd - typo修正、練習問題追加、引用形式変更
- [x] 3-Chapter_2.qmd - チェック完了
- [x] 4-Chapter_3.qmd - チェック完了
- [x] 5-Chapter_4.qmd - チェック完了（set.seed使用済み）
- [x] 6-Chapter_5.qmd - チェック完了（set.seed使用済み）
- [x] 7-Chapter_6.qmd - チェック完了
- [x] 8-degrees_of_freedom.qmd - チェック完了
- [x] 9-Note.qmd - typo修正
- [x] 10-R.qmd - チェック完了
- [x] 11-Replication.qmd - set.seed追加
- [x] 13-Crosstabs.qmd - チェック完了

### Stanによるベイズ推定（14-Bayes.qmd）
- [x] 独立モデル（Independence Model）
- [x] 一様連関モデル（Uniform Association Model）
- [x] RC(1)モデル
- [x] 3元RC(1)モデル（表5.4のデータ）
- [x] 3φモデル（Model 8相当、gnmでは推定不可）
- [x] Unidiffモデル（Xie 1992）
- [x] gnmとの結果比較
- [x] LOO-CV・WAICによるモデル比較
- [x] 技術ノート作成（notes/stan_gnm_*.md）

## 未公開・保留

### 14-Bayes.qmd
- 状態: _quarto.ymlからコメントアウト、.gitignoreに追加
- 理由: 追加コンテンツとして未公開

## 今後の候補タスク

- [ ] 15-Chapter_5_bayes.qmd の作成検討
- [x] Appendix: LEMについて（16-LEM.qmd 作成済み）
- [ ] Stanモデルの追加テスト
- [ ] ドキュメントの最終レビュー

## 技術メモ

### LEM実行方法（macOS）
```bash
wine LEM95.EXE input.INP output.out
# 入力ファイルはCRLF改行が必要
sed -i '' $'s/$/\r/' filename.INP
```

### Stanの注意点
- 対角セルの完全飽和（Unidiffモデル）
- G²計算の補正式（Σμ≠ΣY問題）
- 符号の不定性（RC(1)モデル）
