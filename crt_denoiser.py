import numpy as np
import tensorflow as tf
import os
import random
import pickle
from tqdm import tqdm


current_file = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file)


# Load the model from a file
model = tf.keras.models.load_model(f'{current_directory}/model/model1_2')
model2 = tf.keras.models.load_model(f'{current_directory}/model/model1_2')


def dna_sequence_to_numeric(dna_sequence_list):
    base_mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    num_sequences = len(dna_sequence_list)
    max_sequence_length = max(len(seq) for seq in dna_sequence_list)

    # Initialize a 2D NumPy array filled with zeros
    numeric_array = np.zeros((num_sequences, max_sequence_length), dtype=int)

    for i, seq in enumerate(dna_sequence_list):
        for j, base in enumerate(seq):
            numeric_array[i, j] = base_mapping.get(base, -1)

    return numeric_array


def denoise_image(image: np.ndarray):
    denoised_image = model.predict(image, verbose=0)
    denoised_image2 = model2.predict(np.round(np.array([denoised_image[0]]) * 3) / 3, verbose=0)

    reshaped_test_data = np.round(denoised_image[0] * 3).reshape(100,178)
    reshaped_test_data2 = np.round(denoised_image2[0] * 3).reshape(100,178)

    return reshaped_test_data, reshaped_test_data2


def int_to_char(alignments: np.ndarray) -> str:
    """Maps the char integers to ACGT"""
    array: np.ndarray = np.vectorize({0: "A", 1: "C", 2: "G", 3: "T"}.get)(alignments)
    return array


def denoiser_main(dataset: list):
    """
    Assume data is a 3D array which corresponds to:
    Dimension 1:    Different files
    Dimension 2:    clone_type, CIGAR_string, seq
    Dimension 3:    DNA Sequence
    """

    # Format data to 2D
    concat = []  # O(n) in space might not be desirable
    for data in dataset:
        concat += data

    # mix it 
    concat = random.sample(concat, len(concat))

    # group data into group of 100 and denoise images
    grouped_list = [concat[i:i + 100] for i in range(0, len(concat), 100)]
    
    denoised_data_set1 = []
    denoised_data_set2 = []

    for data in tqdm(grouped_list):
        if len(data) < 100:  # Need to pad this and make this an option
            break
        
        image = []
        clones = []
        cigars = []

        for sub_data in data:
            clone, _, seq = sub_data
            image.append(seq)
            clones.append(clone)
            cigars.append('')
        
        image_array = np.array([dna_sequence_to_numeric(image)]).astype('float32') / 3
        image_array = image_array.reshape(image_array.shape[0], 100, 178, 1)
        denoised_image1, denoised_image2  = denoise_image(image_array)
        
        denoised_data_set1 += zip(clones, cigars, int_to_char(denoised_image1).tolist())
        denoised_data_set2 += zip(clones, cigars, int_to_char(denoised_image2).tolist())

    # Put data back into each corresponding files
    denoised_data_output1 = dict()
    for clone, _, seq_list in list(denoised_data_set1):
        if clone in denoised_data_output1:
            denoised_data_output1[clone].append([clone, '', ''.join(seq_list)])
        else:
            denoised_data_output1[clone] = [[clone, '', ''.join(seq_list)]]
    
    denoised_data_output2 = dict()
    for clone, _, seq_list in list(denoised_data_set2):
        if clone in denoised_data_output2:
            denoised_data_output2[clone].append([clone, '', ''.join(seq_list)])
        else:
            denoised_data_output2[clone] = [[clone, '', ''.join(seq_list)]]

    # Write to files
    for key in denoised_data_output1:
        file_out = f"denoised_{key}1.pkl".replace(":","")
        with open(file_out, "wb") as file:
            pickle.dump(denoised_data_output1[key], file)\

    for key in denoised_data_output2:
        file_out = f"denoised_{key}2.pkl".replace(":","")
        with open(file_out, "wb") as file:
            pickle.dump(denoised_data_output2[key], file)


if __name__ == "__main__":
    files = ["trimmed_data_DD2.pkl", "trimmed_data_3D7.pkl", "trimmed_data_7G8.pkl"]
    # files = ["trimmed_data_DD2.pkl"]
    dataset = []
    for file in files:
        with open(file, 'rb') as f:
            data = pickle.load(f)
            dataset.append(data)
            

    denoiser_main(dataset)