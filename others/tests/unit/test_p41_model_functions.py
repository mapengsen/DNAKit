"""Real upstream P41 model comparisons; no pretrained model/network required.

Version note: validates local Kipoi models/loaders/writers and genomic-benchmarks
CNN training, evaluation and prediction on three deterministic inputs.
"""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path

import pytest

from dnakit.predictions import kipoi_p41_functions as kp
from dnakit.predictions import genomic_benchmarks_p41_functions as gb
from dnakit.exceptions import ConfigurationError
from dnakit.predictions import evo2_p41_functions as evo


@pytest.fixture
def local_kipoi(tmp_path):
    torch = pytest.importorskip('torch')
    pytest.importorskip('kipoi')
    path = tmp_path / 'local_model'
    path.mkdir()
    (path / 'model.py').write_text('''import torch
from torch import nn
class Tiny(nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = nn.Linear(4, 2)
    def forward(self, x):
        return self.linear(x)
''')
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(23)
        torch.save(torch.nn.Linear(4, 2).state_dict(), path / 'linear.pth')
        torch.save({'linear.' + k: v for k,v in torch.load(path / 'linear.pth', weights_only=True).items()}, path / 'weights.pth')
    (path / 'model.yaml').write_text('''type: pytorch
args:
  weights: weights.pth
  module_class: Tiny
  module_file: model.py
  auto_use_cuda: false
info:
  authors:
    - name: DNAKit local test
  doc: Fixed local linear model for API comparison
  license: MIT
schema:
  inputs:
    shape: (4,)
    doc: DNA base composition
  targets:
    shape: (2,)
    doc: Local test outputs
''')
    (path / 'dataloader.py').write_text('''import numpy as np
from kipoi.data import Dataset
class LocalDataset(Dataset):
    def __init__(self, n=5): self.n=n
    def __len__(self): return self.n
    def __getitem__(self, i):
        return {'inputs': np.eye(4, dtype=np.float32)[i % 4], 'targets': np.array([i % 2], dtype=np.float32)}
''')
    (path / 'dataloader.yaml').write_text('''type: Dataset
defined_as: dataloader.py::LocalDataset
args:
  n:
    doc: Number of local records
    optional: true
info:
  authors:
    - name: DNAKit local test
  doc: Local deterministic DNA composition rows
output_schema:
  inputs:
    shape: (4,)
    doc: DNA base composition
  targets:
    shape: (1,)
    doc: Test labels
''')
    return path


@pytest.mark.parametrize('row', [[1,0,0,0], [0.25,0.25,0.25,0.25], [0,0,0,0]])
def test_kipoi_activation_and_gradient_match_real_upstream(local_kipoi, row):
    np = pytest.importorskip('numpy')
    kipoi = pytest.importorskip('kipoi')
    data = np.array([row], dtype=np.float32)
    direct = kipoi.get_model(str(local_kipoi), source='dir', with_dataloader=False)
    expected_activation = direct.predict_activation_on_batch(data, 'linear')
    expected_gradient = direct.input_grad(data, filter_idx=0, layer='linear')
    actual_activation = kp.kipoi_layer_activations(str(local_kipoi), data, 'linear', verbose=False)
    actual_gradient = kp.kipoi_input_gradients(str(local_kipoi), data, filter_idx=0, layer='linear', verbose=False)
    np.testing.assert_allclose(actual_activation['activations'], expected_activation, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(actual_gradient['gradients'], expected_gradient, rtol=1e-6, atol=1e-9)
    json.dumps(actual_activation, allow_nan=False)
    json.dumps(actual_gradient, allow_nan=False)


def test_kipoi_json_float_input_is_accepted(local_kipoi):
    result = kp.kipoi_layer_activations(str(local_kipoi), [[0.25,0.25,0.25,0.25]], 'linear', verbose=False)
    assert result['status'] == 'completed'


@pytest.mark.parametrize('mode', ['raw', 'train', 'predict'])
@pytest.mark.parametrize('n', [1,4,5])
def test_kipoi_loader_modes_and_truncation_match_upstream(local_kipoi, mode, n):
    kipoi = pytest.importorskip('kipoi')
    direct = kipoi.get_dataloader_factory(str(local_kipoi), source='dir')(n=n)
    method = {'raw': direct.batch_iter, 'train': direct.batch_train_iter, 'predict': direct.batch_predict_iter}[mode]
    expected = list(method(batch_size=2, **({'cycle':False} if mode=='train' else {})))
    actual = kp.kipoi_dataloader_batches(str(local_kipoi), dataloader_args={'n':n}, mode=mode, batch_size=2, max_batches=2, verbose=False)
    assert actual['batches'] == kp._plain(expected[:2])
    assert actual['truncated'] == (len(expected)>2)
    assert actual['batch_count'] == min(len(expected),2)


@pytest.mark.parametrize('n', [1,2,3])
def test_kipoi_bed_predictions_match_upstream(tmp_path, n):
    np = pytest.importorskip('numpy')
    pytest.importorskip('kipoi')
    from kipoi import specs, writers
    ranges = {'chr':['chr1']*n, 'start':list(range(n)), 'end':list(range(1,n+1)), 'id':['r'+str(i) for i in range(n)], 'strand':['+']*n}
    preds = np.arange(n*2).reshape(n,2) / 10
    output = tmp_path/'actual.bed'
    direct = tmp_path/'direct.bed'
    schema={'ranges':specs.MetadataStruct(doc='Intervals',type=specs.MetadataType.GENOMIC_RANGES)}
    writer=writers.BedBatchWriter(str(direct),schema)
    writer.batch_write({'metadata':{'ranges':{k:np.asarray(v) for k,v in ranges.items()}},'preds':preds})
    writer.close()
    actual=kp.kipoi_write_bed_predictions(ranges,preds.tolist(),str(output),verbose=False)
    assert output.read_text()==direct.read_text()
    assert len(output.read_text().splitlines()[0].split('\t')) == 8
    assert actual['outputs']==[str(output)]


@pytest.mark.parametrize('values', [[0.0], [0.1,0.2,0.3], [-1.0,0.0,1.0]])
@pytest.mark.parametrize('kind', ['bedgraph','bigwig'])
def test_kipoi_region_track_matches_upstream(tmp_path, values, kind):
    np=pytest.importorskip('numpy')
    pytest.importorskip('kipoi')
    from kipoi import writers
    regions=[{'chr':'chr1','start':0,'end':len(values)}]
    actual=tmp_path/('actual.'+kind)
    direct=tmp_path/('direct.'+kind)
    kwargs={'chrom_sizes':{'chr1':20}} if kind=='bigwig' else {}
    writer=(writers.BedGraphWriter(str(direct)) if kind=='bedgraph' else writers.BigWigWriter(str(direct),chrom_sizes=[('chr1',20)]))
    writer.region_write(regions[0],np.array(values))
    writer.close()
    kp.kipoi_write_region_track(regions,[values],str(actual),output_format=kind,verbose=False,**kwargs)
    if kind=='bedgraph':
        assert actual.read_text()==direct.read_text()
    else:
        pyBigWig=pytest.importorskip('pyBigWig')
        with pyBigWig.open(str(actual)) as a, pyBigWig.open(str(direct)) as e:
            np.testing.assert_allclose(a.values('chr1',0,len(values)),e.values('chr1',0,len(values)),rtol=1e-6,atol=1e-9)


@pytest.mark.parametrize('kind', ['mixed','homopolymer','multiclass'])
def test_genomic_cnn_train_evaluate_predict_match_upstream(tmp_path, kind):
    torch=pytest.importorskip('torch')
    np=pytest.importorskip('numpy')
    pytest.importorskip('genomic_benchmarks')
    from genomic_benchmarks.models.torch import CNN
    torch.set_num_threads(1)
    cfg={'number_of_classes':3 if kind=='multiclass' else 2,'vocab_size':5,'embedding_dim':2,'input_len':100}
    x=np.array([[0]*100,[1]*100,[2]*100] if kind=='homopolymer' else [[(i+j)%5 for j in range(100)] for i in range(3)])
    labels=[0,1,2] if kind=='multiclass' else [0,1,0]
    y=torch.tensor(labels,dtype=torch.float32).reshape(-1,1)
    loader=torch.utils.data.DataLoader(torch.utils.data.TensorDataset(torch.tensor(x,dtype=torch.long),y),batch_size=3,shuffle=False)
    checkpoint=tmp_path/'cnn.pth'
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(71)
        direct=CNN(**cfg,device='cpu')
        with contextlib.redirect_stdout(io.StringIO()): direct.fit(loader,1)
    result=gb.train_genomic_benchmark_cnn(x.tolist(),labels,cfg,str(checkpoint),seed=71,batch_size=3,verbose=False)
    actual_state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    for key,value in direct.state_dict().items():
        torch.testing.assert_close(actual_state[key],value,rtol=1e-6,atol=1e-9)
    direct.eval()
    with torch.inference_mode(): expected_scores=direct(torch.tensor(x,dtype=torch.long)).tolist()
    scores=gb.predict_genomic_benchmark_cnn(x.tolist(),cfg,str(checkpoint),verbose=False)
    np.testing.assert_allclose(scores['scores'],expected_scores,rtol=1e-6,atol=1e-9)
    log=io.StringIO()
    with contextlib.redirect_stdout(log): accuracy,f1=direct.test(loader)
    evaluation=gb.evaluate_genomic_benchmark_cnn(x.tolist(),labels,cfg,str(checkpoint),batch_size=3,verbose=False)
    assert evaluation['accuracy']==accuracy
    assert evaluation['f1']==f1
    assert evaluation['evaluation_log']==log.getvalue()
    json.dumps([result,scores,evaluation],allow_nan=False)


def test_model_rejects_invalid_track_and_cnn_data(tmp_path):
    with pytest.raises(ConfigurationError):
        kp.kipoi_write_region_track([{'chr':'chr1','start':0,'end':1}],[[float('nan')]],str(tmp_path/'bad.bedgraph'),verbose=False)
    with pytest.raises(ConfigurationError):
        gb.train_genomic_benchmark_cnn([[0]], [0], {}, str(tmp_path/'bad.pth'),verbose=False)


@pytest.fixture
def evo_contract(monkeypatch):
    np = pytest.importorskip('numpy')
    torch = pytest.importorskip('torch')
    calls = []
    class Model:
        tokenizer = object()
        model = object()
        def score_sequences(self, sequences, **kwargs):
            calls.append(('score', sequences, kwargs))
            return [-float(len(s)) for s in sequences]
        def forward(self, ids, **kwargs):
            calls.append(('forward', ids.tolist(), kwargs))
            return (torch.zeros((1, 2, 4)), None), {'layer':torch.ones((1, 2, 2))} if kwargs['return_embeddings'] else None
    class Scoring:
        @staticmethod
        def prepare_batch(sequences, tokenizer, **kwargs):
            calls.append(('prepare', sequences, kwargs))
            return torch.tensor([[1, 2]]), [len(s) for s in sequences]
        @staticmethod
        def positional_entropies(sequences, model, tokenizer, **kwargs):
            calls.append(('entropy', sequences, kwargs))
            return [np.arange(len(s), dtype=float) / 10 for s in sequences]
    monkeypatch.setattr(evo, '_load', lambda *args: Model())
    monkeypatch.setattr(evo, 'require', lambda name, extra: np if name=='numpy' else Scoring)
    monkeypatch.setattr(evo, 'version', lambda name: 'contract-only-no-real-weights')
    return calls


@pytest.mark.parametrize('sequence', ['AC','acgt','NNNN'])
def test_evo2_score_contract_preserves_case_and_flags(evo_contract, sequence):
    actual=evo.score_evo2_sequences([sequence], 'unused.ckpt', average_reverse_complement=True, reduce_method='sum', prepend_bos=True, verbose=False)
    assert actual['scores']==[-float(len(sequence))]
    assert evo_contract[-1]==('score',[sequence],{'batch_size':1,'prepend_bos':True,'reduce_method':'sum','average_reverse_complement':True})


@pytest.mark.parametrize('sequence', ['AC','acgt','NNNN'])
def test_evo2_forward_contract_handles_cache_and_embeddings(evo_contract, sequence):
    actual=evo.forward_evo2_sequences([sequence], 'unused.ckpt', layer_names=['layer'], device='cpu', verbose=False)
    assert actual['sequence_lengths']==[len(sequence)]
    assert actual['logits']==[[[0.0]*4]*2]
    assert actual['embeddings']=={'layer':[[[1.0]*2]*2]}
    json.dumps(actual,allow_nan=False)


@pytest.mark.parametrize('sequence', ['AC','acgt','NNNN'])
@pytest.mark.parametrize('metric', ['entropy','perplexity'])
def test_evo2_uncertainty_contract_strand_average(evo_contract, sequence, metric):
    np=pytest.importorskip('numpy')
    actual=evo.evo2_positional_uncertainty(sequence,'unused.ckpt',metric=metric,reverse_complement=True,device='cpu',verbose=False)
    expected=np.full(len(sequence),(len(sequence)-1)/20)
    if metric=='perplexity': expected=np.exp(expected)
    np.testing.assert_allclose(actual['values'],expected,rtol=1e-6,atol=1e-9)
    rc=sequence.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
    assert evo_contract[-1][1]==[sequence,rc]


def test_evo2_invalid_input_and_output_limit_fail_before_model_load(monkeypatch):
    def forbidden(*args): raise AssertionError('model must not load')
    monkeypatch.setattr(evo,'_load',forbidden)
    for sequences in ([], ['A'], ['AC-X'], 'ACGT'):
        with pytest.raises(ConfigurationError): evo.score_evo2_sequences(sequences,'unused',verbose=False)
    with pytest.raises(ConfigurationError):
        evo.forward_evo2_sequences(['AC'],'unused',max_output_values=1,verbose=False)


def test_evo2_missing_local_checkpoint_never_downloads(tmp_path):
    with pytest.raises(Exception, match='does not exist|not found|exist'):
        evo.score_evo2_sequences(['AC'],str(tmp_path/'absent.ckpt'),verbose=False)
