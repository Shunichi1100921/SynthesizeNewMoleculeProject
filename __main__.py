import molecule_synthesizer.controller.synthesis_multi_fragment


if __name__ == '__main__':
    # fragment = input('What do you want to synthesize?').split()
    fragment = ['F1', 'F3']
    molecule_synthesizer.controller.synthesis_multi_fragment.synthesize_fragment()
