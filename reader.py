import time
import serial
from ast import literal_eval
# print("hello")
ser = serial.Serial('/dev/serial0', baudrate=9600,
                    parity=serial.PARITY_NONE, bytesize=serial.EIGHTBITS, timeout=1)


codes = {
    'battery_state': b'\xDD\xA5\x04\x00\xFF\xFD\x77',
    'per_cell_voltage': b'\xDD\xA5\x04\x00\xFF\xFC\x77',
    'bms_version': b'\xDD\xA5\x04\x00\xFF\xFC\x77'
}



def get_cell_data(code):
    print('inside get cell')
    # print('fire:', codes[code])
    try:
        ser.write(codes[code])
        data = ser.readline()
        print(data)
    except:
        print('invalid code')

    
    # print(data[4],data[5])
    #volt1 = data[4] << 8
    #volt1 = volt1|data[5]
    # print(volt1)
    #volt2 = data[6] << 8
    #volt2 = volt2|data[7]
    # print(volt2)
    #print(data, len(data))
    if code == 'per_cell_voltage':
        print('inside per_Cell')
        voltages_ordered = []
        if data:
            print("inside")
            for a in range(4, len(data), 2):
                # print(data[a])
                if (a+1) >= len(data)-2:
                    break
                int_data = (data[a] << 8) | data[a+1]
                float_data = int_data/1000
                print(float_data)
                voltages_ordered.append(float_data)
        print("----------------")

while 1:
    get_cell_data('per_cell_voltage')