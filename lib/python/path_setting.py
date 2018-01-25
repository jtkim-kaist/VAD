
class PathSetting(object):

    def __init__(self, prj_dir, model):
        save_dir = prj_dir + '/data/feat'
        train_dir = save_dir + '/train'
        valid_dir = save_dir + '/valid'
        logs_dir = prj_dir + '/logs/' + model

        self.logs_dir = logs_dir
        self.initial_logs_dir = logs_dir
        self.input_dir = train_dir
        self.output_dir = train_dir+'/Labels'
        self.norm_dir = train_dir
        self.valid_file_dir = valid_dir

