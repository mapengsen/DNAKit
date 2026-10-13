"""Local candidate adapters for the ml models capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def attribute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'attribute', input, options=options)

def build(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'build', input, options=options)

def calculate_prediction_uncertainty(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'calculate_prediction_uncertainty', input, options=options)

def calibrate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'calibrate', input, options=options)

def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'check', input, options=options)

def encode(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'encode', input, options=options)

def evaluate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'evaluate', input, options=options)

def explain(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'explain', input, options=options)

def export_benchmark(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_benchmark', input, options=options)

def export_calibration_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_calibration_model', input, options=options)

def export_chromatin_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_chromatin_model', input, options=options)

def export_cnn_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_cnn_model', input, options=options)

def export_dataset(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_dataset', input, options=options)

def export_expression_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_expression_model', input, options=options)

def export_graph_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_graph_model', input, options=options)

def export_interpretability(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_interpretability', input, options=options)

def export_kmer_embedding(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_kmer_embedding', input, options=options)

def export_leaderboard(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_leaderboard', input, options=options)

def export_long_context_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_long_context_model', input, options=options)

def export_model_onnx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_model_onnx', input, options=options)

def export_onehot_embedding(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_onehot_embedding', input, options=options)

def export_onnx_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_onnx_model', input, options=options)

def export_pangenome_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_pangenome_model', input, options=options)

def export_phenotype_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_phenotype_model', input, options=options)

def export_protein_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_protein_model', input, options=options)

def export_rnn_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_rnn_model', input, options=options)

def export_splice_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_splice_model', input, options=options)

def export_test_split(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_test_split', input, options=options)

def export_tokenizer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_tokenizer', input, options=options)

def export_torchscript_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_torchscript_model', input, options=options)

def export_train_split(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_train_split', input, options=options)

def export_transformer_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_transformer_model', input, options=options)

def export_uncertainty_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_uncertainty_model', input, options=options)

def export_validation_split(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_validation_split', input, options=options)

def export_variant_effect_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'export_variant_effect_model', input, options=options)

def finetune(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'finetune', input, options=options)

def generate_sequence_embeddings(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'generate_sequence_embeddings', input, options=options)

def load(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'load', input, options=options)

def mask_dna_tokens(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'mask_dna_tokens', input, options=options)

def mutagenesis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'mutagenesis', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'predict', input, options=options)

def pretrain(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'pretrain', input, options=options)

def quantize_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'quantize_model', input, options=options)

def run_cross_species_transfer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'run_cross_species_transfer', input, options=options)

def search(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'search', input, options=options)

def serve(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'serve', input, options=options)

def train(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'train', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('ml_models', 'validate', input, options=options)

__all__ = ['attribute', 'build', 'calculate_prediction_uncertainty', 'calibrate', 'check', 'encode', 'evaluate', 'explain', 'export_benchmark', 'export_calibration_model', 'export_chromatin_model', 'export_cnn_model', 'export_dataset', 'export_expression_model', 'export_graph_model', 'export_interpretability', 'export_kmer_embedding', 'export_leaderboard', 'export_long_context_model', 'export_model_onnx', 'export_onehot_embedding', 'export_onnx_model', 'export_pangenome_model', 'export_phenotype_model', 'export_protein_model', 'export_rnn_model', 'export_splice_model', 'export_test_split', 'export_tokenizer', 'export_torchscript_model', 'export_train_split', 'export_transformer_model', 'export_uncertainty_model', 'export_validation_split', 'export_variant_effect_model', 'finetune', 'generate_sequence_embeddings', 'load', 'mask_dna_tokens', 'mutagenesis', 'predict', 'pretrain', 'quantize_model', 'run_cross_species_transfer', 'search', 'serve', 'train', 'validate']

__candidate_adapter_operations__ = set(['attribute', 'build', 'calculate_prediction_uncertainty', 'calibrate', 'check', 'encode', 'evaluate', 'explain', 'export_benchmark', 'export_calibration_model', 'export_chromatin_model', 'export_cnn_model', 'export_dataset', 'export_expression_model', 'export_graph_model', 'export_interpretability', 'export_kmer_embedding', 'export_leaderboard', 'export_long_context_model', 'export_model_onnx', 'export_onehot_embedding', 'export_onnx_model', 'export_pangenome_model', 'export_phenotype_model', 'export_protein_model', 'export_rnn_model', 'export_splice_model', 'export_test_split', 'export_tokenizer', 'export_torchscript_model', 'export_train_split', 'export_transformer_model', 'export_uncertainty_model', 'export_validation_split', 'export_variant_effect_model', 'finetune', 'generate_sequence_embeddings', 'load', 'mask_dna_tokens', 'mutagenesis', 'predict', 'pretrain', 'quantize_model', 'run_cross_species_transfer', 'search', 'serve', 'train', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
