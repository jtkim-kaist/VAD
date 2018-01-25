import graph_save as gs

prj_dir = '/home/sbie/storage3/github/VAD_Toolkit/VAD'

gs.freeze_graph(prj_dir + '/saved_model/temp/temp_LSTM', prj_dir + '/saved_model/graph/LSTM', 'model_1/soft_pred,model_1/raw_labels')

# gs.freeze_graph(prj_dir + '/saved_model/temp_ACAM', prj_dir + '/saved_model/ACAM', 'model_1/logits,model_1/raw_labels')
