def sort_and_average(array, key_array):
    import numpy as np

    povelikosti = np.argsort(key_array)

    array1 = []
    for i in povelikosti:
        array1.append(array[i])
    array = np.array(array1)
    key_array = np.partition(key_array, povelikosti)

    new_array = []
    new_key_array = []
    counts = []
    last = key_array[0]
    new_array.append(array[0])
    new_key_array.append(key_array[0])
    counts.append(1)
    j = 0
    k = 0
    for i in key_array:
        if i > last and j > 0:
            new_key_array.append(i)
            new_array.append(array[j])
            counts.append(1)
            last = i
            k += 1
        elif j > 0:
            new_array[k] += array[j]
            counts[k] += 1
            last = i
        j += 1

    new_array = np.array(new_array)
    new_key_array = np.array(new_key_array)
    counts = np.array(counts)

    new_array = new_array / counts
    return new_array, new_key_array
