# Utils used with tensorflow implementation

import tensorflow as tf
import numpy as np
import scipy.misc as misc
import os, sys
from six.moves import urllib
import tarfile
import zipfile
import scipy.io
import re
import data_reader_bDNN as dr
import data_reader_DNN as dnn_dr

__author__ = 'Juntae'


def vad_test(m_eval, sess_eval, batch_size_eval, eval_file_dir, norm_dir, data_len, eval_type):

    eval_input_dir = eval_file_dir
    eval_output_dir = eval_file_dir + '/Labels'

    pad_size = batch_size_eval - data_len % batch_size_eval
    if eval_type != 2:
        eval_data_set = dr.DataReader(eval_input_dir, eval_output_dir, norm_dir, w=19, u=9, name="eval", pad=pad_size)
    else:
        eval_data_set = dnn_dr.DataReader(eval_input_dir, eval_output_dir, norm_dir, w=19, u=9, name="eval", pad=pad_size)

    final_softout, final_label = evaluation(m_eval, eval_data_set, sess_eval, batch_size_eval, eval_type)

    return final_softout, final_label


def affine_transform(x, output_dim, seed=0, name=None):
    """
    affine transformation Wx+b
    assumes x.shape = (batch_size, num_features)
    """
    initializer = tf.truncated_normal_initializer(stddev=0.02, seed=seed)

    # weights = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim],
    #                           initializer=tf.contrib.layers.xavier_initializer(seed=seed))
    weights = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim],
                              initializer=initializer)
    b = tf.get_variable(name + "_b", [output_dim], initializer=tf.constant_initializer(0.0))

    return tf.matmul(x, weights) + b


def evaluation(m_valid, valid_data_set, sess, eval_batch_size, eval_type):
    # num_samples = valid_data_set.num_samples
    # num_batches = num_samples / batch_size

    if eval_type == 0:  # proposed
        final_softout = []
        final_label = []

        while True:

            valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

            feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: valid_labels,
                         m_valid.keep_probability: 1}
            if valid_data_set.eof_checker():
                final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
                final_label = np.reshape(np.asarray(final_label), [-1, 1])
                valid_data_set.reader_initialize()
                # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
                break

            valid_soft_result, valid_raw_labels = sess.run([m_valid.soft_result, m_valid.raw_labels],
                                                           feed_dict=feed_dict)
            final_softout.append(valid_soft_result)
            final_label.append(valid_raw_labels)

            # if valid_data_set.eof_checker():
            #     final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
            #     final_label = np.reshape(np.asarray(final_label), [-1, 1])
            #     valid_data_set.reader_initialize()
            #     # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            #     break

        return final_softout, final_label

    elif eval_type == 1:  # bdnn
        final_softout = []
        final_label = []

        while True:

            valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

            feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: valid_labels,
                         m_valid.keep_probability: 1}
            if valid_data_set.eof_checker():
                final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
                final_label = np.reshape(np.asarray(final_label), [-1, 1])
                valid_data_set.reader_initialize()
                # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
                break

            valid_cost, valid_logits = sess.run([m_valid.cost, m_valid.logits], feed_dict=feed_dict)
            valid_pred, soft_pred = bdnn_prediction(eval_batch_size + 2*valid_data_set._w, valid_logits, threshold=0.6)
            # print(np.sum(valid_pred))


            raw_indx = int(np.floor(valid_labels.shape[1] / 2))
            raw_labels = valid_labels[:, raw_indx]
            raw_labels = raw_labels.reshape((-1, 1))


            final_softout.append(soft_pred)
            final_label.append(raw_labels)

            # if valid_data_set.eof_checker():
            #     final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
            #     final_label = np.reshape(np.asarray(final_label), [-1, 1])
            #     valid_data_set.reader_initialize()
            #     # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            #     break

        return final_softout, final_label

    elif eval_type == 2:  # dnn

        final_softout = []
        final_label = []

        while True:

            valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

            one_hot_labels = valid_labels.reshape((-1, 1))
            one_hot_labels = dense_to_one_hot(one_hot_labels, num_classes=2)
            feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: one_hot_labels,
                         m_valid.keep_probability: 1}
            if valid_data_set.eof_checker():
                final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
                final_label = np.reshape(np.asarray(final_label), [-1, 1])
                valid_data_set.reader_initialize()
                # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
                break



            soft_pred, raw_labels = sess.run([m_valid.softpred, m_valid.raw_labels], feed_dict=feed_dict)
            raw_labels = raw_labels.reshape((-1, 1))


            final_softout.append(soft_pred)
            final_label.append(raw_labels)

            # if valid_data_set.eof_checker():
            #     final_softout = np.reshape(np.asarray(final_softout), [-1, 1])
            #     final_label = np.reshape(np.asarray(final_label), [-1, 1])
            #     valid_data_set.reader_initialize()
            #     # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            #     break

        return final_softout, final_label



def onehot_tensor(label_batch, num_labels):
    batch_size = label_batch.get_shape().as_list()[0]
    num_labels = tf.cast(num_labels, tf.int32)
    sparse_labels = tf.cast(tf.reshape(label_batch, [-1, 1]), dtype=tf.int32)
    indices = tf.reshape(tf.range(0, batch_size, 1), [batch_size, 1])
    concated = tf.concat(axis=1, values=[indices, sparse_labels])
    outshape = tf.stack([batch_size, num_labels])
    labels = tf.sparse_to_dense(concated, outshape, 1.0, 0.0)
    return labels


def get_model_data(dir_path, model_url):
    maybe_download_and_extract(dir_path, model_url)
    filename = model_url.split("/")[-1]
    filepath = os.path.join(dir_path, filename)
    if not os.path.exists(filepath):
        raise IOError("VGG Model not found!")
    data = scipy.io.loadmat(filepath)
    return data


def maybe_download_and_extract(dir_path, url_name, is_tarfile=False, is_zipfile=False):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    filename = url_name.split('/')[-1]
    filepath = os.path.join(dir_path, filename)
    if not os.path.exists(filepath):
        def _progress(count, block_size, total_size):
            sys.stdout.write(
                '\r>> Downloading %s %.1f%%' % (filename, float(count * block_size) / float(total_size) * 100.0))
            sys.stdout.flush()

        filepath, _ = urllib.request.urlretrieve(url_name, filepath, reporthook=_progress)
        print()
        statinfo = os.stat(filepath)
        print('Succesfully downloaded', filename, statinfo.st_size, 'bytes.')
        if is_tarfile:
            tarfile.open(filepath, 'r:gz').extractall(dir_path)
        elif is_zipfile:
            with zipfile.ZipFile(filepath) as zf:
                zip_dir = zf.namelist()[0]
                zf.extractall(dir_path)


def save_image(image, save_dir, name, mean=None):
    """
    Save image by unprocessing if mean given else just save
    :param mean:
    :param image:
    :param save_dir:
    :param name:
    :return:
    """
    if mean:
        image = unprocess_image(image, mean)
    misc.imsave(os.path.join(save_dir, name + ".png"), image)


def get_variable(weights, name):
    init = tf.constant_initializer(weights, dtype=tf.float32)
    var = tf.get_variable(name=name, initializer=init,  shape=weights.shape)
    return var


def weight_variable(shape, stddev=0.02, name=None):
    # print(shape)
    initial = tf.truncated_normal(shape, stddev=stddev)
    #initial = tf.contrib.layers.xavier_initializer_conv2d()
    if name is None:
        return tf.Variable(initial)
    else:
        return tf.get_variable(name, initializer=initial)


def bias_variable(shape, name=None):
    initial = tf.constant(0.0, shape=shape)
    if name is None:
        return tf.Variable(initial)
    else:
        return tf.get_variable(name, initializer=initial)


def get_tensor_size(tensor):
    from operator import mul
    return reduce(mul, (d.value for d in tensor.get_shape()), 1)


def conv2d_basic(x, W, bias, stride=1):
    conv = tf.nn.conv2d(x, W, strides=[1, stride, stride, 1], padding="SAME")
    return tf.nn.bias_add(conv, bias)


def conv2d_basic_VALID(x, W, bias, stride=1):
    conv = tf.nn.conv2d(x, W, strides=[1, stride, stride, 1], padding="VALID")
    return tf.nn.bias_add(conv, bias)


def conv2d_strided(x, W, b):
    conv = tf.nn.conv2d(x, W, strides=[1, 2, 2, 1], padding="SAME")
    return tf.nn.bias_add(conv, b)


def conv2d_transpose_strided(x, W, b, output_shape=None, stride = 2):
    # print x.get_shape()
    # print W.get_shape()
    if output_shape is None:
        output_shape = x.get_shape().as_list()
        output_shape[1] *= 2
        output_shape[2] *= 2
        output_shape[3] = W.get_shape().as_list()[2]
    # print output_shape
    conv = tf.nn.conv2d_transpose(x, W, output_shape, strides=[1, stride, stride, 1], padding="SAME")
    return tf.nn.bias_add(conv, b)


def leaky_relu(x, alpha=0.0, name=""):
    return tf.maximum(alpha * x, x, name)


def max_pool_2x2(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding="SAME")


def max_pool_2x1(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 1, 1], strides=[1, 2, 1, 1], padding="SAME")


def avg_pool_2x2(x):
    return tf.nn.avg_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding="SAME")


def local_response_norm(x):
    return tf.nn.lrn(x, depth_radius=5, bias=2, alpha=1e-4, beta=0.75)


def batch_norm(x, n_out, phase_train, scope='bn', decay=0.9, eps=1e-5):
    """
    Code taken from http://stackoverflow.com/a/34634291/2267819
    """
    with tf.variable_scope(scope):
        beta = tf.get_variable(name='beta', shape=[n_out], initializer=tf.constant_initializer(0.0)
                               , trainable=True)
        gamma = tf.get_variable(name='gamma', shape=[n_out], initializer=tf.random_normal_initializer(1.0, 0.02),
                                trainable=True)
        batch_mean, batch_var = tf.nn.moments(x, [0, 1, 2], name='moments')
        ema = tf.train.ExponentialMovingAverage(decay=decay)

        def mean_var_with_update():
            ema_apply_op = ema.apply([batch_mean, batch_var])
            with tf.control_dependencies([ema_apply_op]):
                return tf.identity(batch_mean), tf.identity(batch_var)

        mean, var = tf.cond(phase_train,
                            mean_var_with_update,
                            lambda: (ema.average(batch_mean), ema.average(batch_var)))
        normed = tf.nn.batch_normalization(x, mean, var, beta, gamma, eps)
    return normed


def process_image(image, mean_pixel):
    return image - mean_pixel


def unprocess_image(image, mean_pixel):
    return image + mean_pixel


def bottleneck_unit(x, out_chan1, out_chan2, down_stride=False, up_stride=False, name=None):
    """
    Modified implementation from github ry?!
    """

    def conv_transpose(tensor, out_channel, shape, strides, name=None):
        out_shape = tensor.get_shape().as_list()
        in_channel = out_shape[-1]
        kernel = weight_variable([shape, shape, out_channel, in_channel], name=name)
        shape[-1] = out_channel
        return tf.nn.conv2d_transpose(x, kernel, output_shape=out_shape, strides=[1, strides, strides, 1],
                                      padding='SAME', name='conv_transpose')

    def conv(tensor, out_chans, shape, strides, name=None):
        in_channel = tensor.get_shape().as_list()[-1]
        kernel = weight_variable([shape, shape, in_channel, out_chans], name=name)
        return tf.nn.conv2d(x, kernel, strides=[1, strides, strides, 1], padding='SAME', name='conv')

    def bn(tensor, name=None):
        """
        :param tensor: 4D tensor input
        :param name: name of the operation
        :return: local response normalized tensor - not using batch normalization :(
        """
        return tf.nn.lrn(tensor, depth_radius=5, bias=2, alpha=1e-4, beta=0.75, name=name)

    in_chans = x.get_shape().as_list()[3]

    if down_stride or up_stride:
        first_stride = 2
    else:
        first_stride = 1

    with tf.variable_scope('res%s' % name):
        if in_chans == out_chan2:
            b1 = x
        else:
            with tf.variable_scope('branch1'):
                if up_stride:
                    b1 = conv_transpose(x, out_chans=out_chan2, shape=1, strides=first_stride,
                                        name='res%s_branch1' % name)
                else:
                    b1 = conv(x, out_chans=out_chan2, shape=1, strides=first_stride, name='res%s_branch1' % name)
                b1 = bn(b1, 'bn%s_branch1' % name, 'scale%s_branch1' % name)

        with tf.variable_scope('branch2a'):
            if up_stride:
                b2 = conv_transpose(x, out_chans=out_chan1, shape=1, strides=first_stride, name='res%s_branch2a' % name)
            else:
                b2 = conv(x, out_chans=out_chan1, shape=1, strides=first_stride, name='res%s_branch2a' % name)
            b2 = bn(b2, 'bn%s_branch2a' % name, 'scale%s_branch2a' % name)
            b2 = tf.nn.relu(b2, name='relu')

        with tf.variable_scope('branch2b'):
            b2 = conv(b2, out_chans=out_chan1, shape=3, strides=1, name='res%s_branch2b' % name)
            b2 = bn(b2, 'bn%s_branch2b' % name, 'scale%s_branch2b' % name)
            b2 = tf.nn.relu(b2, name='relu')

        with tf.variable_scope('branch2c'):
            b2 = conv(b2, out_chans=out_chan2, shape=1, strides=1, name='res%s_branch2c' % name)
            b2 = bn(b2, 'bn%s_branch2c' % name, 'scale%s_branch2c' % name)

        x = b1 + b2
        return tf.nn.relu(x, name='relu')


def add_to_regularization_and_summary(var):
    if var is not None:
        tf.summary.histogram(var.op.name, var)
        tf.add_to_collection("reg_loss", tf.nn.l2_loss(var))


def add_activation_summary(var):
    if var is not None:
        tf.summary.histogram(var.op.name + "/activation", var)
        tf.summary.scalar(var.op.name + "/sparsity", tf.nn.zero_fraction(var))


def add_gradient_summary(grad, var):
    if grad is not None:
        tf.summary.histogram(var.op.name + "/gradient", grad)


def get_conv_shape(name):
    spec = re.split(':|, |->', name)
    kernel_size = int(spec[5])
    stride = int(spec[7])
    input_fm = int(spec[9])
    output_fm = int(spec[10])
    conv_shape = [kernel_size, kernel_size, input_fm, output_fm]
    return conv_shape, stride


def get_1d_conv_shape(name):
    spec = re.split(':|, |->', name)
    kernel_size = int(spec[5])
    stride = int(spec[7])
    input_fm = int(spec[9])
    output_fm = int(spec[10])
    conv_shape = [kernel_size, 1, input_fm, output_fm]
    return conv_shape, stride


def write_val_summary(graph, loss):
    with graph.as_default():
        val_loss = tf.placeholder(tf.float32, shape=[1], name="loss")
        tf.summary.scalar("entropy", val_loss)
        summary_op = tf.summary.merge_all()
        return summary_op


def conv2lstm_layer(inputs, num_fm):
    """
    make the conv_out flat for rnn input
    :param inputs:
    :param num_fm: # final output feature maps.
    :return: outputs: flattened output. shape = (batch_size, num_fm)
    """
    shape = inputs.get_shape().as_list()
    W = weight_variable([shape[1], shape[2], shape[3], num_fm], name="last_conv_w")
    b = bias_variable([num_fm], name="last_conv_b")
    conv_last = conv2d_basic_VALID(inputs, W, b)
    outputs = tf.nn.relu(conv_last, name="last_relu")
    return outputs


def batch_norm_affine_transform(x, output_dim, decay=0, name=None, seed=0, is_training=True):
    """
    affine transformation Wx+b
    assumes x.shape = (batch_size, num_features)
    """
    # initializer = tf.contrib.layers.xavier_initializer(seed=seed)

    w = tf.get_variable(name+"_w", [x.get_shape()[1], output_dim], initializer = tf.contrib.layers.xavier_initializer(seed=seed))
    b = tf.get_variable(name+"_b", [output_dim], initializer=tf.constant_initializer(0.0))
    affine_result = tf.matmul(x, w) + b
    batch_norm_result = tf.contrib.layers.batch_norm(affine_result, decay=decay, is_training=is_training,
                                                     updates_collections=None)
    return batch_norm_result


def bdnn_transform(inputs, w, u):

    # """
    # :param inputs. shape = (batch_size, feature_size)
    # :param w : decide neighbors
    # :param u : decide neighbors
    # :return: trans_inputs. shape = (batch_size, feature_size*len(neighbors))
    # """

    neighbors_1 = np.arange(-w, -u, u)
    neighbors_2 = np.array([-1, 0, 1])
    neighbors_3 = np.arange(1+u, w+1, u)

    neighbors = np.concatenate((neighbors_1, neighbors_2, neighbors_3), axis=0)

    pad_size = 2*w + inputs.shape[0]
    pad_inputs = np.zeros((pad_size, inputs.shape[1]))
    pad_inputs[0:inputs.shape[0], :] = inputs

    trans_inputs = [np.roll(pad_inputs, -1*neighbors[i], axis=0)[0:inputs.shape[0], :]
                    for i in range(neighbors.shape[0])]

    trans_inputs = np.asarray(trans_inputs)
    trans_inputs = np.transpose(trans_inputs, [1, 0, 2])
    trans_inputs = np.reshape(trans_inputs, (trans_inputs.shape[0], -1))

    return trans_inputs


def bdnn_prediction(bdnn_batch_size, logits, threshold=0.6, w=19, u=9):

    result = np.zeros((bdnn_batch_size, 1))
    indx = np.arange(bdnn_batch_size) + 1
    indx = indx.reshape((bdnn_batch_size, 1))
    indx = bdnn_transform(indx, w, u)
    indx = indx[w:(bdnn_batch_size-w), :]
    indx_list = np.arange(w, bdnn_batch_size - w)

    for i in indx_list:
        indx_temp = np.where((indx-1) == i)
        pred = logits[indx_temp]
        pred = np.sum(pred)/pred.shape[0]
        result[i] = pred + np.random.rand(1)*1e-4

    result = np.trim_zeros(result)
    soft_result = np.float32(result)
    result = result >= threshold

    return result.astype(int), soft_result


def clipped_relu(x, name=None):
    b = tf.get_variable(name+'proposed', [1], initializer=tf.constant_initializer(-.5))
    x = tf.maximum((x+0.5), 0)
    x = tf.minimum(x, 1)

    return x


def dense_to_one_hot(labels_dense, num_classes=2):
    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[(index_offset + labels_dense.ravel()).astype(int)] = 1
    return labels_one_hot.astype(np.float32)