# Contributing to Clinical Genomics Pipeline

We welcome contributions to the Clinical Genomics Pipeline! This document provides guidelines for contributing.

## 🤝 How to Contribute

### Reporting Issues
- Use GitHub Issues to report bugs or request features
- Provide detailed information including:
  - Pipeline version
  - Nextflow version
  - Container engine (Docker/Apptainer)
  - Error messages and logs
  - Sample data (if applicable)

### Submitting Changes
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Run the test suite (`./tests/run_tests.sh`)
6. Commit with descriptive messages
7. Push to your fork
8. Submit a Pull Request

## 🧪 Testing Guidelines

### Unit Tests
```bash
# Run all tests
./tests/run_tests.sh

# Run specific test categories  
python -m pytest tests/test_modules.py::TestPipelineStructure -v
```

### Integration Tests
```bash
# Test with synthetic data
./scripts/run_pipeline.sh --profile test

# Test with stub run (no actual execution)
nextflow run main.nf -profile test -stub-run
```

### Adding New Tests
- Add unit tests for new modules in `tests/test_modules.py`
- Test both success and failure scenarios
- Include edge cases and boundary conditions

## 📋 Development Standards

### Code Style
- Follow Nextflow DSL2 best practices
- Use descriptive process and variable names
- Include comprehensive documentation
- Add version tracking for all tools

### Module Structure
```
modules/
├── category/
│   ├── main.nf          # Workflow definition
│   ├── process1.nf      # Individual processes  
│   └── process2.nf
```

### Container Requirements
- Use open-source tools with permissive licenses
- Build from source or official repositories
- Include version pinning
- Document all dependencies

### Documentation
- Update README.md for major changes
- Add inline comments for complex logic
- Update configuration documentation
- Include usage examples

## 🚀 Release Process

### Version Numbering
We use Semantic Versioning (SemVer):
- **Major** (X.0.0): Breaking changes
- **Minor** (0.X.0): New features, backwards compatible
- **Patch** (0.0.X): Bug fixes, backwards compatible

### Release Checklist
- [ ] All tests passing
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] Version bumped in manifest
- [ ] Tagged release created

## 📞 Getting Help

- **GitHub Discussions**: General questions and usage help
- **GitHub Issues**: Bug reports and feature requests
- **Documentation**: Check `docs/` directory
- **Examples**: See `data/test/` for sample usage

## 🙏 Recognition

Contributors will be acknowledged in:
- README.md contributors section
- Release notes
- Git commit history

Thank you for helping improve clinical genomics workflows! 🧬