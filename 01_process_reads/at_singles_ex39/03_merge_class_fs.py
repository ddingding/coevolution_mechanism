# concatenate classified read files for 2 machines

from constants import (
    AT_SINGLE_CALLED_DIR_M1,
    AT_SINGLE_CALLED_DIR_M2,
    AT_SINGLE_CALLED_DIR_BOTH,
)
from os import listdir

from pipelineTools import concat_files


at_fs1 = sorted(
    [f for f in listdir(AT_SINGLE_CALLED_DIR_M1) if f.endswith("class.csv")]
)
at_fs2 = sorted(
    [f for f in listdir(AT_SINGLE_CALLED_DIR_M2) if f.endswith("class.csv")]
)

list_fs_concat = list(zip(at_fs1, at_fs2))

for fs in list_fs_concat:
    f1, f2 = fs
    f_id = f1[: len("180716Lau_D18-6083")]
    # print(f1, f2)
    assert f1[: len("180716Lau_D18-6083")] == f2[: len("180716Lau_D18-6083")]

    concat_files(
        [AT_SINGLE_CALLED_DIR_M1 + f1, AT_SINGLE_CALLED_DIR_M2 + f2],
        AT_SINGLE_CALLED_DIR_BOTH + f_id + "_class.csv",
    )
