""" Full assembly of the parts to form the complete network """

import torch.nn.functional as F

from .unet_parts import *

class UNet(nn.Module):
    def __init__(self, n_channels, n_classes, bilinear=True):
        super(UNet, self).__init__()
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.bilinear = bilinear

        self.inc = DoubleConv(n_channels, 64)
        self.down1 = Down(64, 128)
        self.down2 = Down(128, 256)
        self.down3 = Down(256, 512)
        factor = 2 if bilinear else 1
        self.down4 = Down(512, 1024 // factor)
        self.up1 = Up(1024, 512 // factor, bilinear)
        self.up2 = Up(512, 256 // factor, bilinear)
        self.up3 = Up(256, 128 // factor, bilinear)
        self.up4 = Up(128, 64, bilinear)
        self.outc = OutConv(64, n_classes)

    def forward(self, x):
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)
        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        x = self.up3(x, x2)
        x = self.up4(x, x1)
        logits = self.outc(x)
        return logits

# Assumes: 1 channel image, 1 class
class MyUNet(nn.Module):
    def __init__(self, n_size=32):
        super(MyUNet, self).__init__()
        bilinear = True
        n_channels = 1  # Grayscale images
        n_classes = 1   # We have pupil (255) + background (0)

        self.inc = DoubleConv(n_channels, n_size)
        self.down1 = Down(n_size, n_size*2)
        self.down2 = Down(n_size*2, n_size*4)
        self.down3 = Down(n_size*4, n_size*8)
        self.down4 = Down(n_size*8, n_size*8)

        self.up1 = Up(n_size*16, n_size*4, bilinear)
        self.up2 = Up(n_size*8, n_size*2, bilinear)
        self.up3 = Up(n_size*4, n_size, bilinear)
        self.up4 = Up(n_size*2, n_size, bilinear)
        self.outc = OutConv(n_size, n_classes)

    def forward(self, x):
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)
        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        x = self.up3(x, x2)
        x = self.up4(x, x1)
        logits = self.outc(x)
        return logits
