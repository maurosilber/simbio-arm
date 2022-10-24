# simbio-arm

## Installation

```
pip install simbio-arm
```

## Usage

```python
from simbio.models.arm import ARM_extrinsic
from simbio.simulator import Simulator

t = range(100)
Simulator(ARM_extrinsic).run(t).plot()
```
