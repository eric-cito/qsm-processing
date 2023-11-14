from qsm.QSMPipeline import QSMPipeline
from qsm.FileLocations import FileLocations

locs = FileLocations('/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/'); # must end with /s)
pipeline = QSMPipeline(locs)
pipeline.Run()