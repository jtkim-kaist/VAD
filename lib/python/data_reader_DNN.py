import numpy as np
import os
import glob
import utils
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


class DataReader(object):

    def __init__(self, input_dir, output_dir, norm_dir, w=19, u=9, name=None, pad=None):
        # print(name.title() + " data reader initialization...")
        self._input_dir = input_dir
        self._output_dir = output_dir
        self._norm_dir = norm_dir
        self._input_file_list = sorted(glob.glob(input_dir+'/*.bin'))
        self._input_spec_list = sorted(glob.glob(input_dir+'/*.txt'))
        self._output_file_list = sorted(glob.glob(output_dir+'/*.bin'))
        self._file_len = len(self._input_file_list)
        self._name = name
        assert self._file_len == len(self._output_file_list), "# input files and output file is not matched"

        self._num_file = 0
        self._start_idx = 0
        if pad is not None:
            self._inputs = self._read_input(self._input_file_list[self._num_file], self._input_spec_list[self._num_file])
            self._inputs = np.concatenate((self._inputs, np.zeros((pad, self._inputs.shape[1]), dtype=np.float32)))
            self._outputs = self._read_output(self._output_file_list[self._num_file])
            self._outputs = np.concatenate((self._outputs, np.zeros((pad, self._outputs.shape[1]), dtype=np.float32)))
        else:
            self._inputs = self._read_input(self._input_file_list[self._num_file], self._input_spec_list[self._num_file])
            self._outputs = self._read_output(self._output_file_list[self._num_file])
        self._w = w
        self._u = u
        self.eof = False
        self.file_change = False
        self._outputs = self._outputs[0:self._inputs.shape[0]]
        self._num_figure = 1
        assert np.shape(self._inputs)[0] == np.shape(self._outputs)[0], \
            ("# samples is not matched between input: %d and output: %d files"
             % (np.shape(self._inputs)[0], np.shape(self._outputs)[0]))

        self.num_samples = np.shape(self._outputs)[0]

        norm_param = sio.loadmat(self._norm_dir+'/global_normalize_factor.mat')
        self.train_mean = norm_param['global_mean']
        self.train_std = norm_param['global_std']

        # print("Done")
        # print("BOF : " + self._name + " file_" + str(self._num_file).zfill(2))

    def _binary_read_with_shape(self):
        pass

    @staticmethod
    def _read_input(input_file_dir, input_spec_dir):

        data = np.fromfile(input_file_dir, dtype=np.float32)  # (# total frame, feature_size)
        with open(input_spec_dir,'r') as f:
            spec = f.readline()
            size = spec.split(',')
        data = data.reshape((int(size[0]), int(size[1])), order='F')

        return data

    @staticmethod
    def _read_output(output_file_dir):

        data = np.fromfile(output_file_dir, dtype=np.float32)  # data shape : (# total frame,)
        data = data.reshape(-1, 1)  # data shape : (# total frame, 1)

        return data

    def next_batch(self, batch_size):

        if self._start_idx + batch_size > self.num_samples:

            self._start_idx = 0
            self.file_change = True
            self._num_file += 1

            # print("EOF : " + self._name + " file_" + str(self._num_file-1).zfill(2) +
            #       " -> BOF : " + self._name + " file_" + str(self._num_file).zfill(2))

            if self._num_file > self._file_len - 1:
                self.eof = True
                self._num_file = 0
                # print("EOF : last " + self._name + " file. " + "-> BOF : " + self._name + " file_" +
                #       str(self._num_file).zfill(2))

            self._inputs = self._read_input(self._input_file_list[self._num_file], self._input_spec_list[self._num_file])
            self._outputs = self._read_output(self._output_file_list[self._num_file])

            data_len = np.shape(self._inputs)[0]
            self._outputs = self._outputs[0:data_len, :]

            assert np.shape(self._inputs)[0] == np.shape(self._outputs)[0], \
                ("# samples is not matched between input: %d and output: %d files"
                 % (np.shape(self._inputs)[0], np.shape(self._outputs)[0]))

            self.num_samples = np.shape(self._outputs)[0]
            # print("current file number : %d, samples : %d" % (self._num_file + 1, self.num_samples))
            #print("Loaded " + self._name + " file number : %d" % (self._num_file + 1))
        else:
            self.file_change = False
            self.eof = False

        inputs = self._inputs[self._start_idx:self._start_idx + batch_size, :]
        inputs = self.normalize(inputs)
        inputs = utils.bdnn_transform(inputs, self._w, self._u)
        inputs = inputs[self._w: (batch_size-self._w), :]

        outputs = self._outputs[self._start_idx:self._start_idx + batch_size, :]

        # if valid:
        #     plt.figure(self._num_figure)
        #     self._num_figure += 1
        #     bb = np.zeros(aa.shape)
        #     bb[:, 200:250] = outputs*10
        #
        #     cc = aa + bb
        #     imgplot = plt.imshow(cc.T)
        #     plt.show()

        outputs = outputs[self._w: (batch_size - self._w), :]

        self._start_idx += batch_size
        # print(self._start_idx)
        # print(self.num_samples)
        return inputs, outputs

        #num_batches = (np.shape(self._outputs)[0] - np.shape(self._outputs)[0] % batch_size) / batch_size
    def normalize(self, x):
        x = (x - self.train_mean)/self.train_std
        # a = (np.std(x, axis=0))
        return x

    def reader_initialize(self):
        self._num_file = 0
        self._start_idx = 0
        self.eof = False

    def eof_checker(self):
        return self.eof

    def file_change_checker(self):
        return self.file_change

    def file_change_initialize(self):
        self.file_change = False


def dense_to_one_hot(labels_dense, num_classes=2):
    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[index_offset + labels_dense.ravel()] = 1
    return labels_one_hot


# file_dir = "/home/sbie/github/VAD_KJT/Datamake/Database/Aurora2withSE"
# input_dir1 = file_dir + "/STFT2"
# output_dir1 = file_dir + "/Labels"
# dr = DataReader(input_dir1, output_dir1, input_dir1,name='test')
#
# for i in range(1000000):
#     tt, pp = dr.next_batch(500)
#     print("asdf")



