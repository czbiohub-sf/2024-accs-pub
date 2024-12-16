import logging


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import accs24_pub_figs.scripts as scripts
    scripts.generate_all_figs()
