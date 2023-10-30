from typing import Annotated, Literal

import numpy as np
import numpy.typing as npt

LArray3 = Annotated[npt.NDArray[np.float64], Literal[3]]
LArray3_3 = Annotated[npt.NDArray[np.float64], Literal[3,3]]
LArrayN = Annotated[npt.NDArray[np.float64], Literal['N']]
LArrayN_3 = Annotated[npt.NDArray[np.float64], Literal['N',3]]
