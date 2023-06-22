#!/usr/bin/env python

import sys
import io
import struct
from dataclasses import dataclass
import array
import itertools
from enum import Enum, auto

PEAK_SENSITIVITY = 0.001
HISTOGRAM_CLOCK_TOLERANCE = 0.05

STARTING_CLOCK_WEIGHT = 0.6
CURRENT_CLOCK_WEIGHT = 0.3
PULSE_CLOCK_WEIGHT = 1.0 - (STARTING_CLOCK_WEIGHT + CURRENT_CLOCK_WEIGHT)
PULSE_CLOCK_TOLERANCE = 0.1
STARTING_CLOCK_NO_PULSE_WEIGHT = 0.65
CURRENT_CLOCK_NO_PULSE_WEIGHT = 0.35

class RecordFormat(Enum):
    MFM = auto()

class RangeMethod(Enum):
    LARGEST_GAP = auto(),
    MIDDLE_LOWEST = auto(),
    MIDDLE = auto()

@dataclass(frozen=True)
class TrackIndex:
    pos : int
    sample_counter : int
    index_time : int

graph_bar = array.array('u', itertools.repeat('#', 70))

@dataclass
class HistogramInfo:
    histogram : (None, {int: int}) = None
    peaks : (None, [int]) = None
    troughs : (None, [[int]]) = None
    ranges : (None, [int], [(int, int)]) = None
    peak_threshold : int = 0.0

    def print_histogram(self):
        keys = sorted(self.keys())
        highest = 0
        for key in keys:
            if self[key] > highest:
                highest = self[key]
        for key in keys:
            print("{:03} {:05} ".format(key, self[key]), end='')
            barlen = int(self[key] / (highest / len(graph_bar)))
            if barlen > 0:
                print(graph_bar[:barlen-1].tounicode())
            else:
                print()

@dataclass
class ClockInfo:
    starting_clock : float = 0.0
    record_format : RecordFormat = None
    stream : 'StreamFile' = None
    current_clock : float = 0.0
    index : int = -1
    buffer_pos : int = 0
    midCycle : bool = False
    lastShort : int = 0

    def reset(self):
        self.midCycle = False
        self.lastShort = 0

    def goto_index(self, index):
        self.buffer_pos = self.stream.indexes[index].pos
        self.index = index

    def interpret_pulse(self, reverse, last_time = 0.0, last_pulse = 0.0):
        this_time = last_time
        this_pulse = last_pulse + self.stream.buffer[self.buffer_pos]
        clocks = 0
        while this_time < this_pulse - (self.current_clock * 0.5):
            this_time += self.current_clock
            clocks += 1
        if clocks <= 1 and self.record_format == RecordFormat.MFM:
            # MFM can't have a pulse 1 clock long so try to string
            # multiple together assuming there was an extra pulse
            if last_time == 0.0:
                print("WARNING: 1 clock transition in MFM fmroat.")
            if reverse:
                if self.buffer_pos == 0:
                    self.buffer_pos = len(self.stream.buffer)
                    self.index = len(self.stream.indexes) - 1
                self.buffer_pos -= 1
            else:
                self.buffer_pos += 1
                if self.buffer_pos == len(self.stream.buffer):
                    self.buffer_pos = 0
                    self.index = -1
            clocks = self.interpret_pulse(reverse, this_time, this_pulse)
        newclock = self.stream.buffer[self.buffer_pos] / clocks
        if newclock > self.current_clock * (1.0 - PULSE_CLOCK_TOLERANCE) and \
           newclock < self.current_clock * (1.0 + PULSE_CLOCK_TOLERANCE):
            self.current_clock = ((self.starting_clock * STARTING_CLOCK_WEIGHT) +
                                  (self.current_clock * CURRENT_CLOCK_WEIGHT) +
                                  (newclock * PULSE_CLOCK_WEIGHT))
        else:
            self.current_clock = ((self.starting_clock * STARTING_CLOCK_NO_PULSE_WEIGHT) +
                                  (self.current_clock * CURRENT_CLOCK_NO_PULSE_WEIGHT))
        return clocks

    def interpret_pulses(self, num=-1, reverse=False, until_end=False, until_index=False):
        buffer = array.array('B', itertools.repeat(0, num))
        if num < 0:
            until_end = True
        count = 0
        while num < 0 or count < num:
            if reverse:
                if self.buffer_pos == 0:
                    self.buffer_pos = len(self.stream.buffer)
                    self.index = len(self.stream.indexes) - 1
                    if until_index or until_end:
                        return buffer
                self.buffer_pos -= 1
                if self.index >= 0 and \
                   self.buffer_pos < self.stream.indexes[self.index-1].pos:
                    self.index -= 1
                    if until_index:
                        return buffer
            buffer.append(self.interpret_pulse(reverse))
            if not reverse:
                self.buffer_pos += 1
                if self.index < len(self.stream.indexes) - 1 and \
                   self.buffer_pos > self.stream.indexes[self.index+1].pos:
                    self.index += 1
                    if until_index or until_end:
                        return buffer
                if self.buffer_pos == len(self.stream.buffer):
                    self.buffer_pos = 0
                    self.index = -1
                    if until_index:
                        return buffer
            count += 1
        return buffer

    def interpret_mfm(self, inBuffer, num = -1):
        outBuffer = array.array('B')
        pos = 0
        count = 0
        while (num < 0 or pos < num) and count < len(inBuffer):
            # == 1 should be impossible, but it'l still be caught, just with a rather odd message
            if inBuffer[count] == 2:
                if self.midCycle:
                    outBuffer.append(1)
                else:
                    outBuffer.append(0)
                pos += 1
            elif inBuffer[count] == 3:
                if self.midCycle:
                    outBuffer.append(1)
                    outBuffer.append(0)
                    pos += 2
                    self.midCycle = False
                else:
                    outBuffer.append(0)
                    pos += 1
                    self.midCycle = True
            elif inBuffer[count] == 4:
                if self.midCycle:
                    outBuffer.append(1)
                    outBuffer.append(0)
                else:
                    print("WARNING: 4 clock pulse where uninterpretable? Maybe assume a pair of 2 clock pulses? {} @ {:X}({:X})".format(inBuffer[count], self.buffer_pos, count))
                    outBuffer.append(0)
                    outBuffer.append(0)
                pos += 2
            else:
                print("WARNING: Overly long pulse? {} > {} @ {:X}".format(inBuffer[count], self.current_clock * 4.0, self.buffer_pos))
                if inBuffer[count] % 2 == 1:
                    self.midCycle = not self.midCycle
            count += 1
        return outBuffer

def bsearch(buffer, val):
    high = len(buffer)-1
    low = 0
    lasti = -1
    while True:
        i = low + ((high - low) // 2)
        if i == lasti:
            raise IndexError("Failed to find value {}.".format(val))
        lasti = i
        if val < buffer[i]:
            high = i
        elif val > buffer[i]:
            low = i
        else:
            return i

class StreamFile:
    def __init__(self, filename, peak_sensitivity=PEAK_SENSITIVITY):
        self.buffer = array.array('I')
        with open(filename, 'rb') as infile:
            self.peak_sensitivity = peak_sensitivity
            self.sck = None
            self.ick = None
            self.indexes = []
            # Read in all out of band data and collect info
            last_index = None
            offsets = array.array('I')
            self.totalval = 0
            val = 0
            while True:
                blockhdr, = struct.unpack('B', infile.read(1))
                if blockhdr <= 0x7:
                    offsets.append(infile.tell()-1)
                    val += blockhdr * 0x100 + infile.read(1)[0]
                    self.buffer.append(val)
                    self.totalval += val
                    val = 0
                elif blockhdr == 0x08:
                    pass
                elif blockhdr == 0x09:
                    infile.read(1)
                elif blockhdr == 0x0a:
                    infile.read(2)
                elif blockhdr == 0x0b:
                    val += 0x10000
                elif blockhdr == 0x0c:
                    offsets.append(infile.tell()-1)
                    val += infile.read(1)[0] * 0x100 + infile.read(1)[0]
                    self.buffer.append(val)
                    self.totalval += val
                    val = 0
                elif blockhdr == 0x0d:
                    oob_type, = struct.unpack('B', infile.read(1))
                    oob_size, = struct.unpack('H', infile.read(2))
                    match oob_type:
                        case 0x00:
                            print("Invalid Block")
                        case 0x01:
                            stream_pos, transfer_time = struct.unpack('II', infile.read(oob_size))
                            print("Stream Info: Pos: {}, Time: {}".format(stream_pos, transfer_time))
                        case 0x02:
                            stream_pos, sample_counter, index_counter = struct.unpack('III', infile.read(oob_size))
                            print("Index: Pos: {}, Sample Counter: {}, Index Counter: {}".format(stream_pos, sample_counter, index_counter))
                            if last_index == None:
                                last_index = index_counter
                                self.indexes.append(TrackIndex(bsearch(offsets, stream_pos),
                                                               sample_counter, 0))
                            else:
                                if index_counter < last_index:
                                    last_index -= 2**32
                                index_diff = index_counter - last_index
                                last_index = index_counter
                                last_time = index_diff/self.ick
                                last_RPM = 60/last_time
                                print("Last index was {} index clocks ago ({} seconds, {} RPM).".format(index_diff, last_time, last_RPM))
                                self.indexes.append(TrackIndex(bsearch(offsets, stream_pos),
                                                               sample_counter,
                                                               self.indexes[-1].index_time+index_diff))
                        case 0x03:
                            stream_pos, status_code = struct.unpack_from('II', infile.read(oob_size))
                            print("Stream End: Pos: {}, Status: {} ".format(stream_pos, status_code), end='')
                            match status_code:
                                case 0x00:
                                    print("OK")
                                case 0x01:
                                    print("Buffer Error")
                                case 0x02:
                                    print("No Index")
                                case unknown:
                                    print("Unknown {}".format(unknown))
                        case 0x04:
                            print("Device Info:")
                            info = infile.read(oob_size)[:-1].decode('cp437')
                            for kv in info.split(', '):
                                k, v = kv.split('=')
                                if k == 'sck':
                                    self.sck = float(v)
                                elif k == 'ick':
                                    self.ick = float(v)
                                print("{}={}".format(k, v))
                        case 0x0d:
                            print("EOF")
                            break
                        case unknown:
                            print("Unknown {}".format(unknown))
                elif blockhdr >= 0x0e:
                    offsets.append(infile.tell()-1)
                    val += blockhdr
                    self.buffer.append(val)
                    self.totalval += val
                    val = 0
        self.histogramInfo, self.clockInfo = self.build_histogram(0, len(self.buffer))

    def build_histogram(self, start, size,
                        range_method=RangeMethod.LARGEST_GAP,
                        histogramInfo=None, clockInfo=None,
                        peak_sensitivity=None):
        if peak_sensitivity is None:
            peak_sensitivity = self.peak_sensitivity
        histogram = {}
        peaks = []
        troughs = []
        ranges = []
        if histogramInfo is not None and histogramInfo.histogram is not None:
            histogram = histogramInfo.histogram
        else:
            val = 0
            bufpos = start
            while bufpos < start + size:
                try:
                    histogram[self.buffer[bufpos]] += 1
                except KeyError:
                    histogram[self.buffer[bufpos]] = 1
                bufpos += 1

        keys = sorted(histogram.keys())

        if histogramInfo is not None and histogramInfo.peak_threshold > 0:
            peak_threshold = histogramInfo.peak_threshold
        else:
            peak_threshold = int(self.totalval / len(keys) * peak_sensitivity)

        if histogramInfo is None or (histogramInfo is not None and
                                     histogramInfo.peaks is not None and
                                     histogramInfo.troughs is not None):
            # find peaks
            minval = keys[0]
            minkey = 0
            maxval = keys[0]
            maxkey = 0
            find_peak = False
            find_trough = False
            for key in keys:
                val = histogram[key]
                if val < minval:
                    minval = val
                    minkey = key
                if val > maxval:
                    maxval = val
                    maxkey = key
                if not find_trough:
                    if not find_peak:
                        if maxval - minval > peak_threshold:
                            troughs.append(minkey)
                            find_peak = True
                    else:
                        if maxval - val > peak_threshold:
                            peaks.append(maxkey)
                            minval = maxval
                            minkey = maxkey
                            find_peak = False
                            find_trough = True
                else:
                    if val - minval > peak_threshold:
                        maxval = minval
                        maxkey = minkey
                        find_trough = False
            # remove entries before peaks
            if troughs[0] < peaks[0]:
                troughs = troughs[1:]
            if troughs[-1] > peaks[-1]:
                troughs = troughs[:-1]
            for i in range(len(troughs)):
                if troughs[i] < peaks[i] or peaks[i+1] < troughs[i]:
                    return histogram, peaks, troughs, None

        if histogramInfo is not None and histogramInfo.peaks is not None:
            peaks = histogramInfo.peaks

        if histogramInfo is not None and histogramInfo.troughs is not None:
            troughs = histogramInfo.troughs

        if not isinstance(troughs[0], list):
            # find all lowest keys between peaks
            for i in range(len(peaks)-1):
                thistrough = []
                for key in keys[keys.index(peaks[i])+1:keys.index(peaks[i+1])]:
                    if histogram[key] == histogram[troughs[i]]:
                        thistrough.append(key)
                troughs[i] = thistrough

        if histogramInfo is not None and histogramInfo.ranges is not None:
            if not isinstance(histogramInfo.ranges, tuple):
                histogramInfo.ranges[0] = (keys[0], histogramInfo.ranges[0])
                for i in range(len(histogramInfo.ranges)-1):
                    histogramInfo.ranges[i+1] = (histogramInfo.ranges[i], histogramInfo.ranges[i+1])
                histogramInfo.ranges[-1] = (histogramInfo.ranges[-1], keys[-1])
        else:
            # find peak ranges
            if range_method == RangeMethod.LARGEST_GAP:
                for i in range(len(peaks)-1):
                    last_key = peaks[i]
                    largest_gap = last_key
                    largest_gap_size = 0
                    multiple_largest_gap = False
                    for key in keys[keys.index(peaks[i])+1:keys.index(peaks[i+1])]:
                        if key - last_key > largest_gap_size:
                            largest_gap = last_key
                            largest_gap_size = key - last_key
                            multiple_largest_gap = False
                        elif key - last_key == largest_gap_size:
                            multiple_largest_gap = True
                        last_key = key
                    if multiple_largest_gap:
                        print("WARNING: No clear histogram gap found, using the middle lowest value found instead!")
                        if i == 0:
                            ranges.append((keys[0], troughs[0][len(troughs[0])//2]))
                        else:
                            ranges.append((keys[keys.index(ranges[-1][1])+1], troughs[i][len(troughs[i])//2]))
                    else:
                        if i == 0:
                            ranges.append((keys[0], largest_gap))
                        else:
                            ranges.append((keys[keys.index(ranges[-1][1])+1], largest_gap))
            else:
                if range_method == RangeMethod.MIDDLE_LOWEST:
                    ranges.append((keys[0], troughs[0][len(troughs[0])//2]))
                elif range_method == RangeMethod.MIDDLE:
                    firstpeak = keys.index(peaks[0])
                    secondpeak = keys.index(peaks[1])
                    ranges.append((keys[0], keys[firstpeak + ((secondpeak - firstpeak) // 2)]))
                for i in range(len(troughs)-1):
                    if range_method == RangeMethod.MIDDLE_LOWEST:
                        ranges.append((keys[keys.index(ranges[-1][1])+1], troughs[i+1][len(troughs[i+1])//2]))
                    elif range_method == RangeMethod.MIDDLE:
                        firstpeak = keys.index(peaks[i+1])
                        secondpeak = keys.index(peaks[i+2])
                        ranges.append((keys[keys.index(ranges[-1][1])+1], keys[firstpeak + ((secondpeak - firstpeak) // 2)]))
            ranges.append((keys[keys.index(ranges[-1][1])+1], keys[-1]))

        if len(peaks) == 3:
            record_format = RecordFormat.MFM
            if clockInfo is None or \
               (clockInfo is not None and clockInfo.starting_clock <= 0):
                starting_clock = float(peaks[0]) / 2.0
                clock2 = float(peaks[1]) / 3.0
                clock3 = float(peaks[2]) / 4.0
                if clock2 > starting_clock * (1.0 - HISTOGRAM_CLOCK_TOLERANCE) and \
                   clock2 < starting_clock * (1.0 + HISTOGRAM_CLOCK_TOLERANCE):
                    if clock3 > starting_clock * (1.0 - HISTOGRAM_CLOCK_TOLERANCE) and \
                       clock3 < starting_clock * (1.0 + HISTOGRAM_CLOCK_TOLERANCE):
                        starting_clock = (starting_clock + clock2 + clock3) / 3.0
                    else:
                        print("WARNING: Clock 3 is {}% out!".format(HISTOGRAM_CLOCK_TOLERANCE * 100.0))
                        starting_clock = (starting_clock + clock2) / 2.0
                else:
                    print("WARNING: Clock 2 is {}% out!".format(HISTOGRAM_CLOCK_TOLERANCE * 100.0))
                    if clock3 > starting_clock * (1.0 - HISTOGRAM_CLOCK_TOLERANCE) and \
                       clock3 < starting_clock * (1.0 + HISTOGRAM_CLOCK_TOLERANCE):
                        starting_clock = (starting_clock + clock3) / 2.0
                    else:
                        print("WARNING: Clock 3 is {}% out!".format(HISTOGRAM_CLOCK_TOLERANCE * 100.0))
        else:
            raise ValueError("Unknown recording format.")

        if histogramInfo is not None:
            if histogramInfo.histogram is None:
                histogramInfo.histogram = histogram
            if histogramInfo.peaks is None:
                histogramInfo.peaks = peaks
            if histogramInfo.troughs is None:
                histogramInfo.troughs = troughs
            if histogramInfo.ranges is None:
                histogramInfo.ranges = ranges
            if histogramInfo.peak_threshold <= 0.0:
                histogramInfo.peak_threshold = peak_threshold
        else:
            histogramInfo = HistogramInfo(histogram, peaks, troughs, ranges, peak_threshold)

        if clockInfo is not None:
            if clockInfo.starting_clock <= 0.0:
                clockInfo.starting_clock = starting_clock
            if clockInfo.stream is None:
                clockInfo.stream = self
            if clockInfo.record_format is None:
                clockInfo.record_format = record_format
        else:
            clockInfo = ClockInfo(starting_clock, record_format, self)

        clockInfo.current_clock = clockInfo.starting_clock

        return histogramInfo, clockInfo 

def print_buffer(buffer):
    for val, i in enumerate(buffer):
        if i % 64 == 0:
            print("\n{:04X}".format(i), end='')
        if i % 8 == 0:
            print(end=' ')
        print(val, end='')
    print()

def print_buffers_compare(buffer1, buffer2):
    for i, (val1, val2) in enumerate(zip(buffer1, buffer2)):
        if i % 64 == 0:
            print("\n{:04X}".format(i), end='')
        if i % 8 == 0:
            print(end=' ')
        if val1 == val2:
            print('.', end='')
        else:
            print(abs(val1 - val2), end='')
    print()

infile = StreamFile(sys.argv[1])
print(infile.histogramInfo.peaks)
print(infile.histogramInfo.troughs)
print(infile.histogramInfo.ranges)
print(infile.histogramInfo.peak_threshold)
print(infile.indexes)
#print_histogram(infile.histogramInfo.histogram)
infile.clockInfo.goto_index(2)
buffer = infile.clockInfo.interpret_pulses(-1, until_index=True)
infile.clockInfo.goto_index(3)
buffer2 = infile.clockInfo.interpret_pulses(-1, until_index=True)
print(len(buffer))
print(len(buffer2))
#print_buffers_compare(buffer, buffer2)
#buffer = infile.clockInfo.interpret_mfm(buffer)
#print_buffer(buffer)
histogramInfo, clockInfo = infile.build_histogram(infile.indexes[0].pos, infile.indexes[1].pos - infile.indexes[0].pos)
print(histogramInfo.ranges)
histogramInfo, clockInfo = infile.build_histogram(infile.indexes[0].pos, infile.indexes[1].pos - infile.indexes[0].pos, range_method=RangeMethod.MIDDLE_LOWEST)
print(histogramInfo.ranges)
histogramInfo, clockInfo = infile.build_histogram(infile.indexes[0].pos, infile.indexes[1].pos - infile.indexes[0].pos, range_method=RangeMethod.MIDDLE)
#print_histogram(histogramInfo.histogram)
print(histogramInfo.peaks)
print(histogramInfo.troughs)
print(histogramInfo.ranges)
print(infile.histogramInfo.peak_threshold)
