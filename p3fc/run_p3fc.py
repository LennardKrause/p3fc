def main():
    import sys
    import logging
    from PyQt6 import QtWidgets
    import p3fc
    from p3fc.lib.classes import Main_GUI

    # create logger
    logging.basicConfig(level=logging.DEBUG, style='{', format='{message:>20s} @ {funcName}()')
    logging.debug(__name__)

    app = QtWidgets.QApplication(sys.argv)
    w = Main_GUI()
    w.setWindowTitle('Convert PILATUS3 Data to Bruker Format, {} | lkrause@chem.au.dk'.format(p3fc.__version__))
    w.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
