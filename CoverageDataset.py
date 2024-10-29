import torch
from torch.utils.data import Dataset
import pyBigWig as bw
import pandas as pd
import os
import numpy as np
import hashlib 
from src.gena_lm.utils import get_service_token_encodings, concatenate_encodings
import logging
import gc
import pickle
import h5py
import tqdm
import tempfile
import sys
from transformers import AutoTokenizer

class CoverageDataset(Dataset):
    N_SERVICE_TOKENS = 2

    def __init__(self,
                 tokenizer,
                 targets_path: str,
                 tokenized_genome_path: str,
                 key_file_path: str,  # Path to the text file with keys
                 subset: float = 1,
                 seed: int = 42,
                 max_seq_len: int = 512,
                 n_context_tokens: int = 0,
                 n_target_tokens: int = 510,
                 shift_length: int = 510,
                 loglevel: int = logging.WARNING,
                 pybedtools_tempdir: str = "/tmp",
                 hash_path: str = None,
                 force_h5: bool = False,
                 transform_targets=None,
                 ):
        """
        Args:
            targets_path: path to bigWig tab-delimeted metadata file, containing for each bigWig label and filepath
            tokenized_genome_path: path to the hdf5 file containing tokenized genome data
            key_file_path: path to the text file containing keys for hdf5 file
            n_context_tokens: number of tokens for left and right context
            n_target_tokens: number of target tokens
            shift_length: shift between samples (in tokens). shift_length==2*n_context_tokens+n_target_tokens for non-overlapping samples
        """        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(level=loglevel)
        self.logger.info("Initializing dataset")
        self.pybedtools_tempdir = pybedtools_tempdir

        if isinstance(tokenizer, str):
            self.tokenizer = AutoTokenizer.from_pretrained(tokenizer)
        else:
            self.tokenizer = tokenizer
        self.max_seq_len = max_seq_len
        self.n_context_tokens = n_context_tokens
        self.n_target_tokens = n_target_tokens
        self.sample_length = 2 * n_context_tokens + n_target_tokens
        assert self.sample_length + self.N_SERVICE_TOKENS == max_seq_len
        self.shift_length = shift_length
        self.service_token_encodings = get_service_token_encodings(self.tokenizer)
        self.transform_targets = transform_targets

        self.read_bigWig_metadata(targets_path=targets_path)
        self.n_targets = self.n_tracks

        self.tokenized_genome_path = tokenized_genome_path
        self.subset = subset
        self.seed = seed

        if hash_path is None:
            self.hash_path = self.get_hash_path()
        else:
            self.hash_path = hash_path


        # Load keys from the key file
        self.load_keys(key_file_path)

    def load_keys(self, key_file_path):
        """Loads the list of keys from a text file."""
        with open(key_file_path, 'r') as f:
            self.keys = [line.strip() for line in f]
        self.logger.info(f"Loaded {len(self.keys)} keys from {key_file_path}")

    def read_bigWig_metadata(self, targets_path):
        self.logger.info("Reading bigWig file paths")
        self.bigWigPaths = {}
        with open(targets_path) as fin:
            for line in fin:
                line = line.strip()
                if line and not line.startswith("#"):
                    k, v = line.split()
                    assert k not in self.bigWigPaths.keys(), \
                        f"Found repeated bigWig name {k}"
                    self.bigWigPaths[k] = v
        self.n_tracks = len(self.bigWigPaths.keys())
        self.files_opened = False

    def open_files(self): #load_tokenized_genome called after forking
        self.load_tokenized_genome()
        self.bigWigHandlers = {}
        for k, v in self.bigWigPaths.items():
            try:
                self.bigWigHandlers[k] = bw.open(v)
            except Exception as e:
                self.logger.error(f"Error opening bigwig file {v}")
                print(e.__traceback__) 
        self.files_opened = True

    def load_tokenized_genome(self):
        """Loads the tokenized genome data from an hdf5 file."""
        self.tokenized_genome_data = {}
        self.tokenized_genome_data =  h5py.File(self.tokenized_genome_path, 'r')
                
    def get_hash_path(self):
        m = hashlib.blake2b(digest_size=8)
        m.update(str(self.tokenized_genome_path).encode("utf-8"))
        m.update(str(self.sample_length).encode("utf-8"))
        if self.subset < 1.0:
            m.update(str(self.subset).encode("utf-8"))
            m.update(str(self.seed).encode("utf-8"))
        self.hash_suffix = m.hexdigest()
        hash_path = str(self.tokenized_genome_path) + "." + self.hash_suffix
        return hash_path

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, idx):
        if not self.files_opened:
            self.open_files()
        
        # Get the key for the current index
        region_key = self.keys[idx]
        
        chrom, start, end = region_key.split(',')
        
        input_ids = self.tokenized_genome_data[region_key][...]             
        features = {
            "input_ids": input_ids.astype(np.int32).reshape(1, -1),
            "attention_mask": np.ones(len(input_ids), dtype=bool).reshape(1, -1),
            "token_type_ids": np.zeros(len(input_ids), dtype=np.int32).reshape(1, -1),
        }
        
        # Concatenate with CLS and SEP tokens
        features = concatenate_encodings(
            [
                self.service_token_encodings["CLS"],
                features,
                self.service_token_encodings["SEP"],
            ]
        )
        
        # Reducing extra dimension
        features = {k: v[0] for k, v in features.items()}
        
        labels = np.zeros(shape=(self.n_tracks,), dtype=np.float32)
        bins_mask = np.ones(shape=(self.sample_length + self.N_SERVICE_TOKENS,), dtype=bool)
        
        bins_mask[0:self.n_context_tokens + 1] = 0
        bins_mask[-1 - self.n_context_tokens:] = 0
        assert sum(bins_mask) == self.n_target_tokens
        
        for ind, (k, v) in enumerate(self.bigWigHandlers.items()):
            try:
                values = v.values(str(chrom), int(start), int(end), numpy=True)
                labels[ind] = np.sum(values)
            except RuntimeError as e:
                print(idx, chrom, int(start), int(end))
                raise e
        
        if self.transform_targets is not None:
            features["labels"] = self.transform_targets(labels)
        else:
            features["labels"] = labels
        features["bins_mask"] = bins_mask
        
        return features

def worker_init_fn(worker_id):
    worker_info = torch.utils.data.get_worker_info()
    dataset = worker_info.dataset
    dataset.open_files()

class logtransform():
    def __init__(self, pseudocount=0.1):
        self.pseudocount = pseudocount
        self.rounddigits = int(abs(np.log10(pseudocount))) + 2

    def __call__(self, x):
        return np.log(x + self.pseudocount)

    def reverse(self, x):
        return np.round(np.exp(x) - self.pseudocount, self.rounddigits)

class logTransformWithPseudocount():
    def __init__(self, pseudocount):
        self.pseudocount = pseudocount
        self.rounddigits = int(abs(np.log10(pseudocount))) + 2

    def __call__(self, x):
        return np.log(x + self.pseudocount)

    def reverse(self, x):
        return np.round(np.exp(x) - self.pseudocount, self.rounddigits)