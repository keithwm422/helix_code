import socket
s = socket.socket()
address = '127.0.0.1'
port = 1025  # port number is a number, not string
s.bind((address,port))
print s.recv(1024)
s.close
