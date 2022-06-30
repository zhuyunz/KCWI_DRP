from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessObject(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # create mask and flag arrays
        if self.action.args.ccddata.flags is None:
            self.action.args.ccddata.flags = np.zeros(
                self.action.args.ccddata.data.shape, dtype=np.uint8)
        if self.action.args.ccddata.mask is None:
            self.action.args.ccddata.mask = np.zeros(
                self.action.args.ccddata.data.shape, dtype=np.uint8)

        # saturation check
        sat = (self.action.args.ccddata.data >= 65535)
        self.action.args.ccddata.flags[sat] += 1
        self.action.args.ccddata.mask[sat] = 1

        return self.action.args
