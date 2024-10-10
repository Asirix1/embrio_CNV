import os
import torch

torch.cuda.empty_cache()


os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import json
import logging
import os
from pathlib import Path

import torch
from torch.utils.data import DataLoader, DistributedSampler
import transformers
from transformers import AutoConfig, HfArgumentParser
from sklearn.metrics import mean_squared_error
import numpy as np

from lm_experiments_tools import Trainer, TrainerArgs, get_optimizer
from lm_experiments_tools.utils import get_cls_by_name, collect_run_configuration, get_git_diff
import lm_experiments_tools.optimizers as optimizers

from CoverageDataset import CoverageDataset, worker_init_fn

import horovod.torch as hvd

from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf
from hydra.utils import instantiate

import time
timestamp = time.strftime("%Y%m%d-%H%M%S")


def batch_transform_fn(batch):
    return {
        'input_ids': batch['input_ids'],
        'attention_mask': batch['attention_mask'],
        'labels': batch['labels']
    }


def metrics_fn(output):
    metrics = {'loss': output['loss'].detach()}
    return metrics

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

if os.environ.get('CUDA_VISIBLE_DEVICES', None) is None:
    os.environ['CUDA_VISIBLE_DEVICES'] = ','.join([str(i) for i in range(torch.cuda.device_count())])

logger.info(f"CUDA_VISIBLE_DEVICES: {os.environ['CUDA_VISIBLE_DEVICES']}")
logger.info(f"CUDA DEVICE COUNT: {torch.cuda.device_count()}")

hvd.init()
torch.set_num_threads(4)
torch.cuda.set_device(hvd.local_rank())

parser = HfArgumentParser(TrainerArgs)
parser.add_argument('--experiment_config', type=str, help='path to the experiment config')
parser.add_argument('--resume', type=str, default=None, help='resume checkpoint file name')

if __name__ == '__main__':
    args = parser.parse_args()
    experiment_config_path = Path(args.experiment_config).expanduser().absolute()
    with initialize_config_dir(str(experiment_config_path.parents[0])):
        experiment_config = compose(config_name=experiment_config_path.name)

    if "args_params" in experiment_config:
        trainer_kwargs = instantiate(experiment_config["args_params"])
        for k, v in trainer_kwargs.items():
            if hasattr(args, k):
                logger.warning(f"Conflicting setting for option {k}, overwritten by config option {v}")
            args.__setattr__(k, v)

    args.model_path = os.path.join(args.model_path, timestamp)

    if hvd.rank() == 0:
        if args.model_path is None:
            raise ValueError("Model path should not be None")

    logger.info(f"hvd.rank(): {hvd.rank()}, Model path: {args.model_path}")

    if hvd.rank() == 0 and args.model_path is not None:
        model_path = Path(args.model_path)
        if not model_path.exists():
            model_path.mkdir(parents=True)
        args_dict = collect_run_configuration(args)
        json.dump(args_dict, open(model_path / 'config.json', 'w'), indent=4)
        open(model_path / 'git.diff', 'w').write(get_git_diff())

    per_worker_batch_size = args.batch_size * args.gradient_accumulation_steps
    kwargs = {'pin_memory': True, 'num_workers': args.data_n_workers}


    logger.info(f'Preparing training data from: {experiment_config["train_dataset"]["tokenized_genome_path"]}')
    train_dataset = instantiate(experiment_config["train_dataset"])
    train_sampler = DistributedSampler(train_dataset, rank=hvd.rank(), num_replicas=hvd.size(), shuffle=True)
    train_dataloader = DataLoader(train_dataset, batch_size=per_worker_batch_size, sampler=train_sampler, worker_init_fn=worker_init_fn, **kwargs)

    valid_dataloader = None
    if "valid_dataset" in experiment_config:
        logger.info(f'Preparing validation data from: {experiment_config["valid_dataset"]["tokenized_genome_path"]}')
        valid_dataset = instantiate(experiment_config["valid_dataset"])
        valid_sampler = DistributedSampler(valid_dataset, rank=hvd.rank(), num_replicas=hvd.size(), shuffle=False)
        valid_dataloader = DataLoader(valid_dataset, batch_size=per_worker_batch_size, sampler=valid_sampler, worker_init_fn=worker_init_fn, **kwargs)
    model_cfg = AutoConfig.from_pretrained(args.model_cfg)
    model_cfg.num_labels = train_dataset.n_targets

    model_cls = get_cls_by_name(args.model_cls)
    model = model_cls(config=model_cfg)


    optimizer_cls = get_optimizer(args.optimizer)
    optimizer = optimizer_cls(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    trainer = Trainer(args, model, optimizer, train_dataloader, valid_dataloader=valid_dataloader, train_sampler=train_sampler, batch_transform_fn=batch_transform_fn, metrics_fn=metrics_fn)
    trainer.train()
    

    if "valid_dataset" in experiment_config:
        logger.info('Running validation on validation data:')
        trainer.validate(valid_dataloader)
