import torch
from torch import nn
from transformers.modeling_outputs import TokenClassifierOutput

from src.gena_lm.modeling_bert import BertPreTrainedModel, BertModel

class BertForCoveragePrediction(BertPreTrainedModel):

    _keys_to_ignore_on_load_unexpected = [r"pooler"]

    def __init__(self, config, 
                        loss_fct = nn.PoissonNLLLoss(log_input=False, reduction='none'),  
                        num_fc_layers = 0,
                        activation = nn.Softplus()):  
        super().__init__(config)
        self.num_labels = config.num_labels
        self.config = config

        self.bert = BertModel(config, add_pooling_layer=False)
        classifier_dropout = (
            config.classifier_dropout if config.classifier_dropout is not None else config.hidden_dropout_prob
        )
        self.dropout = nn.Dropout(classifier_dropout)

        if num_fc_layers > 0:
            self.fc_layers = nn.ModuleList()
            self.ReLU = nn.LeakyReLU()
            for i in range(num_fc_layers):
                self.fc_layers.append(nn.Linear(config.hidden_size, config.hidden_size))
                self.fc_layers.append(self.ReLU)
                self.fc_layers.append(self.dropout)

        self.classifier = nn.Linear(config.hidden_size, config.num_labels)        
        
        self.activation = activation
        self.loss_fct = loss_fct

        # Initialize weights and apply final processing
        self.post_init()

    def forward(
        self,
        input_ids=None,
        attention_mask=None,
        token_type_ids=None,
        position_ids=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
        target_weights=None,
    ):
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict
        outputs = self.bert(
            input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict
        )
        bins_output = outputs.last_hidden_state[:, 0, :]
        print(bins_output.shape)

 
        if hasattr(self,"fc_layers"):
            for layer in self.fc_layers:
                bins_output = layer(bins_output)
        
        logits = self.classifier(bins_output)
        pred = self.activation(logits)  

        loss = None
        if labels is not None:
            if target_weights is not None:
                loss = self.loss_fct(pred, labels)
                assert len(loss.shape) > 1, f"loss shape == 1, {loss.shape}. Did you use reduction='None'?"
                assert loss.shape == target_weights.shape
                loss = (loss * target_weights) // (target_weights != 0).sum()
                loss = loss.mean()
            else:
                loss = self.loss_fct(pred, labels)
                loss = loss.mean()  
        else:
            raise ValueError("labels are none")

        if not return_dict:
            output = (logits,) + outputs[2:]
            return ((loss,) + output) if loss is not None else output

        return TokenClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )